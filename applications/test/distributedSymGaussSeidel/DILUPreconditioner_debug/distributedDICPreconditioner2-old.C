/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019,2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "distributedDICPreconditioner2.H"
#include "processorLduInterface.H"
#include "cyclicAMILduInterface.H"
//#include "processorColour.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedDICPreconditioner2, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<distributedDICPreconditioner2>
        adddistributedDICPreconditioner2SymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distributedDICPreconditioner2::meshColour
(
    const lduMesh& lm,
    labelList& cellColour
) const
{
    const lduAddressing& addr = lm.lduAddr();
    const label* const __restrict__ uPtr = addr.upperAddr().begin();
    const label* const __restrict__ lPtr = addr.lowerAddr().begin();
    const label* const __restrict__ ownStartPtr = addr.ownerStartAddr().begin();

    const label nCells = addr.size();

    cellColour.setSize(nCells, -1);

    if (nCells == 0)
    {
        return;
    }

    label colouri = -1;
    label celli = 0;
    while (true)
    {
        // Find new seed
        for (; celli<nCells; celli++)
        {
            if (cellColour[celli] == -1)
            {
                break;
            }
        }

        if (celli == -1 || celli == nCells)
        {
            break;
        }

        colouri++;

        // Walk front
        DynamicList<label> front;
        front.append(celli);

        DynamicList<label> newFront;
        while (true)
        {
            newFront.clear();
            for (const label celli : front)
            {
                cellColour[celli] = colouri;

                const label fStart = ownStartPtr[celli];
                const label fEnd = ownStartPtr[celli + 1];

                for (label facei=fStart; facei<fEnd; facei++)
                {
                    const label nbr =
                    (
                        lPtr[facei] == celli
                      ? uPtr[facei]
                      : lPtr[facei]
                    );
                    if (cellColour[nbr] == -1)
                    {
                        newFront.append(nbr);
                    }
                }
            }

            if (newFront.empty())
            {
                break;
            }

            front.transfer(newFront);
        }
    }
}


void Foam::distributedDICPreconditioner2::updateMatrixInterfaces
(
    const bool add,
    const FieldField<Field, solveScalar>& coupleCoeffs,
    const labelList& selectedInterfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt
) const
{
    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();
    const auto& interfaces = solver_.interfaces();

    const label startOfRequests = Pstream::nRequests();

    // From lduMatrix::initMatrixInterfaces

    // Avoid any conflicts with inter-processor comms
    UPstream::msgType() += 1;

    for (const label interfacei : selectedInterfaces)
    {
        interfaces[interfacei].initInterfaceMatrixUpdate
        (
            result,
            add,
            lduAddr,
            interfacei,
            psiif,
            coupleCoeffs[interfacei],
            cmpt,
            UPstream::commsTypes::nonBlocking
        );
    }

    UPstream::waitRequests(startOfRequests);

    // Consume
    for (const label interfacei : selectedInterfaces)
    {
        interfaces[interfacei].updateInterfaceMatrix
        (
            result,
            add,
            lduAddr,
            interfacei,
            psiif,
            coupleCoeffs[interfacei],
            cmpt,
            UPstream::commsTypes::nonBlocking
        );
    }

    UPstream::msgType() -= 1;
}


void Foam::distributedDICPreconditioner2::receive
(
    const labelList& selectedInterfaces
) const
{
    const auto& interfaces = solver_.interfaces();
    const auto& interfaceBouCoeffs = solver_.interfaceBouCoeffs();

    // Start reads
    for (const label inti : selectedInterfaces)
    {
        const auto& intf = interfaces[inti].interface();
        const auto* ppp = isA<const processorLduInterface>(intf);

        auto& recvBuf = recvBufs_[inti];
        recvBuf.setSize(interfaceBouCoeffs[inti].size());

        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            ppp->neighbProcNo(),
            recvBuf.data_bytes(),
            recvBuf.size_bytes(),
            ppp->tag(),
            ppp->comm()
        );
    }
}


void Foam::distributedDICPreconditioner2::send
(
    const labelList& selectedInterfaces,
    const solveScalarField& psiInternal
) const
{
    const auto& interfaces = solver_.interfaces();

    // Start writes
    for (const label inti : selectedInterfaces)
    {
        const auto& intf = interfaces[inti].interface();
        const auto* ppp = isA<const processorLduInterface>(intf);
        const auto& faceCells = intf.faceCells();

        auto& sendBuf = sendBufs_[inti];

        sendBuf.setSize(faceCells.size());
        forAll(faceCells, face)
        {
            sendBuf[face] = psiInternal[faceCells[face]];
        }

        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            ppp->neighbProcNo(),
            sendBuf.cdata_bytes(),
            sendBuf.size_bytes(),
            ppp->tag(),
            ppp->comm()
        );
    }
}


void Foam::distributedDICPreconditioner2::sendGlobal
(
    const labelList& selectedInterfaces,
    solveScalarField& psi,
    const label colouri
) const
{
    const auto& interfaces = solver_.interfaces();

    // Save old data
    FieldField<Field, solveScalar> one(interfaces.size());
    FieldField<Field, solveScalar> old(interfaces.size());
    for (const label inti : selectedInterfaces)
    {
        const auto& intf = interfaces[inti].interface();
        const auto& fc = intf.faceCells();
        one.set(inti, new solveScalarField(fc.size(), 1.0));
        old.set(inti, new solveScalarField(psi, intf.faceCells()));
    }
    updateMatrixInterfaces
    (
        false,              // add to psi
        one,
        selectedInterfaces, //global_,
        psi,                // send data
        psi,                // result
        0                   // cmpt
    );

    colourBufs_.setSize(nColours_);
    if (!colourBufs_.set(colouri))
    {
        colourBufs_.set(colouri, new FieldField<Field, solveScalar>(nColours_));
    }
    auto& colourBuf = colourBufs_[colouri];
    for (const label inti : selectedInterfaces)
    {
        const auto& intf = interfaces[inti].interface();
        const auto& fc = intf.faceCells();

        if (!colourBuf.set(inti))
        {
            colourBuf.set(inti, new solveScalarField(fc.size()));
        }
        auto& cb = colourBuf[inti];
        auto& oldValues = old[inti];

        forAll(cb, face)
        {
            const label cell = fc[face];
            // Store change in boundary value
            cb[face] = psi[cell]-oldValues[face];
            // Restore old value
            std::swap(psi[cell], oldValues[face]);
        }
    }
}


void Foam::distributedDICPreconditioner2::check() const
{
    forAll(coeffs_, inti)
    {
        if (coeffs_.set(inti))
        {
            const auto& coeff = coeffs_[inti];
            for (const auto& c : coeff)
            {
                if (c != solveScalar(0.0))
                {
                    FatalErrorInFunction<< "int:" << inti << exit(FatalError);
                }
            }
        }
    }
}


void Foam::distributedDICPreconditioner2::zeroCoeffs
(
    const labelList& selectedInterfaces
) const
{
    // Reset value to zero
    for (const label inti : selectedInterfaces)
    {
        coeffs_[inti] = Zero;
    }
}


void Foam::distributedDICPreconditioner2::calcReciprocalD
(
    solveScalarField& rD
) const
{
    const auto& interfaces = solver_.interfaces();
    const auto& interfaceBouCoeffs = solver_.interfaceBouCoeffs();
    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();

    // Start swapping remote contributions
    const label startOfRequests = Pstream::nRequests();

    // Start reads (into recvBufs_)
    receive(lowerNbrs_);

    // Start with diagonal
    const scalarField& diag = matrix.diag();
    rD.setSize(diag.size());
    std::copy(diag.begin(), diag.end(), rD.begin());

    Pout<< "** calcReciprocalD :"
        << " diag:" << flatOutput(rD) << endl;


    const label* const __restrict__ uPtr = lduAddr.upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr.lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();


    // Subtract coupled contributions
    {
        // Wait for finish. Received result in recvBufs
        UPstream::waitRequests(startOfRequests);

        for (const label inti : lowerNbrs_)
        {
            const auto& intf = interfaces[inti].interface();
            // TBD: do not use patch faceCells but passed-in addressing?
            const auto& faceCells = intf.faceCells();
            const auto& recvBuf = recvBufs_[inti];
            const auto& bc = interfaceBouCoeffs[inti];

            Pout<< "inti:" << inti
                << " field:" << 1.0/recvBuf
                << " inserting into cells:" << flatOutput(faceCells)
                << " with coeffs:" << flatOutput((bc*bc)())
                << endl;

            forAll(recvBuf, face)
            {
                // Note:interfaceBouCoeffs is -upperPtr
                rD[faceCells[face]] -= bc[face]*bc[face]/recvBuf[face];
            }
        }
        Pout<< "** after proc lowerNbrs_:" << flatOutput(rD) << endl;
    }


    for (label colouri = 0; colouri < nColours_; colouri++)
    {
        Pout<< "** starting colour:" << colouri << endl;

        // Do non-processor boundaries
        // Problem is that neighbour patch should
        // pick up result of completed owner patch (so after internal
        // faces) and not the current result
        //{
        //    check();
        //    for (const label inti : lowerGlobal_)
        //    {
        //        const auto& intf = interfaces[inti].interface();
        //        // TBD: do not use patch faceCells but passed-in addressing?
        //        const auto& faceCells = intf.faceCells();
        //        auto& coeff = coeffs_[inti];
        //        const auto& bc = interfaceBouCoeffs[inti];
        //        forAll(coeff, face)
        //        {
        //            const label cell = faceCells[face];
        //            // Receiving side
        //            if (cellColour_[cell] == colouri)
        //            {
        //                coeff[face] = bc[face]*bc[face];
        //            }
        //        }
        //        Pout<< "For int:" << inti << " setting coefficients to "
        //            << flatOutput(coeff) << endl;
        //    }
        //
        //    // We want to subtract bc[face]*bc[face]/recvBuf so
        //    // - coefficient is bc*bc
        //    // - input field is 1/rD
        //    // - subtract from rD
        //    solveScalarField invRD(1.0/rD);
        //    updateMatrixInterfaces
        //    (
        //        true,           // subtract from rD
        //        coeffs_,
        //        lowerGlobal_,   //global_,
        //        invRD,          // send data
        //        rD,             // result
        //        0               // cmpt
        //    );
        //    // Reset value to zero
        //    zeroCoeffs(lowerGlobal_);
        //    Pout<< "** after lowerGlobal_:" << flatOutput(rD) << endl;
        //}
        {
            check();
            for (const label inti : lowerGlobal_)
            {
                const auto& intf = interfaces[inti].interface();
                // TBD: do not use patch faceCells but passed-in addressing?
                const auto& faceCells = intf.faceCells();
                const auto& recvBuf = colourBufs_[colouri][inti];
                const auto& bc = interfaceBouCoeffs[inti];

                forAll(recvBuf, face)
                {
                    // Note:interfaceBouCoeffs is -upperPtr
                    rD[faceCells[face]] -= bc[face]*bc[face]/recvBuf[face];
                }
            }
            Pout<< "** after lowerGlobal_:" << flatOutput(rD) << endl;
        }


        const label nFaces = matrix.upper().size();
        for (label face=0; face<nFaces; face++)
        {
            const label cell = lPtr[face];
            if (cellColour_[cell] == colouri)
            {
                rD[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rD[lPtr[face]];
            }
        }

        // Store effect of exchanging rD to higher interfaces in colourBufs_
        sendGlobal(higherGlobal_, rD, colouri);

        Pout<< "** finished colour:" << colouri << nl << endl;
    }


    // Start writes of rD (using sendBufs)
    send(higherNbrs_, rD);

    // Calculate the reciprocal of the preconditioned diagonal
    const label nCells = rD.size();

    for (label cell=0; cell<nCells; cell++)
    {
        rD[cell] = 1.0/rD[cell];
    }

    // Wait until all finished. Necessary? Cannot interleave reciprocal
    // calc with comms anyway.
    UPstream::waitRequests(startOfRequests);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedDICPreconditioner2::distributedDICPreconditioner2
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD_(sol.matrix().diag().size())
{
    const lduMesh& mesh = sol.matrix().mesh();
    const auto& interfaces = sol.interfaces();
    const auto& interfaceBouCoeffs = sol.interfaceBouCoeffs();


    // Allocate buffers
    // ~~~~~~~~~~~~~~~~

    coeffs_.setSize(interfaces.size());
    sendBufs_.setSize(interfaces.size());
    recvBufs_.setSize(interfaces.size());
    forAll(interfaceBouCoeffs, inti)
    {
        if (interfaces.set(inti))
        {
            const auto& bc = interfaceBouCoeffs[inti];
            sendBufs_.set(inti, new scalarField(bc.size(), Zero));
            recvBufs_.set(inti, new scalarField(bc.size(), Zero));
        }
    }


    // Determine processor colouring
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //const processorColour& colours = processorColour::New(mesh);
    //const label colouri = colours[Pstream::myProcNo()];
    const label colouri = Pstream::myProcNo();

    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            const auto& intf = interfaces[inti].interface();
            const auto* ppp = isA<const processorLduInterface>(intf);
            if (ppp)
            {
                //const label nbrColour = colours[ppp->neighbProcNo()];
                const label nbrColour = ppp->neighbProcNo();
                if (nbrColour < colouri)
                {
                    lowerNbrs_.append(inti);
                }
                else if (nbrColour > colouri)
                {
                    higherNbrs_.append(inti);
                }
                else
                {
                    WarningInFunction
                        << "weird processorLduInterface"
                        << " from " << ppp->myProcNo()
                        << " to " << ppp->neighbProcNo()
                        << endl;
                }
            }
            else
            {
                global_.append(inti);
                const auto* AMIpp = isA<const cyclicAMILduInterface>(intf);
                if (AMIpp)
                {
                    if (!AMIpp->owner())
                    {
                        lowerGlobal_.append(inti);
                    }
                    else
                    {
                        higherGlobal_.append(inti);
                    }
                    const auto& bc = interfaceBouCoeffs[inti];
                    coeffs_.set(inti, new solveScalarField(bc.size(), Zero));
                }
            }
        }
    }
    DebugVar(global_);
    DebugVar(lowerGlobal_);
    DebugVar(higherGlobal_);


    // Determine local colouring/zoneing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    meshColour(mesh, cellColour_);
    nColours_ = max(cellColour_)+1;
    reduce(nColours_, maxOp<label>(), mesh.comm());

    Pout<< "Determining nColours:" << nColours_
        << " cellColour_:" << flatOutput(cellColour_) << endl;

    calcReciprocalD(rD_);
    DebugVar(rD_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distributedDICPreconditioner2::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA,
    const direction
) const
{
    const auto& interfaces = solver_.interfaces();
    const auto& interfaceBouCoeffs = solver_.interfaceBouCoeffs();
    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();

    solveScalar* __restrict__ wAPtr = wA.begin();
    const solveScalar* __restrict__ rAPtr = rA.begin();
    const solveScalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr = lduAddr.upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr.lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    const label nCells = wA.size();
    const label nFaces = matrix.upper().size();
    const label nFacesM1 = nFaces - 1;


    // Forward sweep
    // ~~~~~~~~~~~~~

    label startOfRequests = Pstream::nRequests();

    // Start reads (into recvBufs)
    receive(lowerNbrs_);

    // Initialise 'internal' cells
    for (label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }

    Pout<< "** after scaling with diag:" << flatOutput(wA) << endl;

    // Do 'halo' contributions from lower numbered procs
    {
        // Wait for finish. Received result in recvBufs
        UPstream::waitRequests(startOfRequests);

        for (const label inti : lowerNbrs_)
        {
            const auto& intf = interfaces[inti].interface();
            // TBD: do not use patch faceCells but passed-in
            // addressing?
            const auto& faceCells = intf.faceCells();
            const auto& recvBuf = recvBufs_[inti];
            const auto& bc = interfaceBouCoeffs[inti];

            forAll(recvBuf, face)
            {
                // Note: interfaceBouCoeffs is -upperPtr
                const label cell = faceCells[face];
                wAPtr[cell] += rDPtr[cell]*bc[face]*recvBuf[face];
            }
        }
        Pout<< "** after lowerNbrs_:" << flatOutput(wA) << endl;


        // Do non-processor boundaries
        for (label colouri = 0; colouri < nColours_; colouri++)
        {
            Pout<< "** starting forward colour:" << colouri << endl;

            check();
            for (const label inti : lowerGlobal_)
            {
                // Note: interfaceBouCoeffs is -upperPtr
                const auto& intf = interfaces[inti].interface();
                auto& coeff = coeffs_[inti];
                const auto& faceCells = intf.faceCells();
                const auto& bc = interfaceBouCoeffs[inti];
                forAll(coeff, face)
                {
                    const label cell = faceCells[face];
                    if (cellColour_[cell] == colouri)
                    {
                        coeff[face] = rDPtr[cell]*bc[face];
                    }
                }
            }

            const solveScalarField sendField(wA);
            updateMatrixInterfaces
            (
                false,          // add to rD
                coeffs_,
                lowerGlobal_,   //global_,
                sendField,      // send data
                wA,             // result
                0               // cmpt
            );
            // Reset value to zero
            zeroCoeffs(lowerGlobal_);
            Pout<< "** after lowerGlobal_:" << flatOutput(wA) << endl;

            for (label face=0; face<nFaces; face++)
            {
                const label cell = lPtr[face];
                if (cellColour_[cell] == colouri)
                {
                    wAPtr[uPtr[face]] -=
                        rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
                }
            }

            Pout<< "** finished forward colour:" << colouri << nl << endl;
        }
    }

    Pout<< "** after forward sweep:" << flatOutput(wA) << endl;

    // Start writes of wA (using sendBufs)
    send(higherNbrs_, wA);


    // Backward sweep
    // ~~~~~~~~~~~~~~

    // Start receives
    receive(higherNbrs_);

    {
        UPstream::waitRequests(startOfRequests);

        for (const label inti : higherNbrs_)
        {
            const auto& intf = interfaces[inti].interface();
            // TBD: do not use patch faceCells but passed-in
            // addressing?
            const auto& faceCells = intf.faceCells();
            const auto& recvBuf = recvBufs_[inti];
            const auto& bc = interfaceBouCoeffs[inti];

            forAll(recvBuf, face)
            {
                // Note: interfaceBouCoeffs is -upperPtr
                const label cell = faceCells[face];
                wAPtr[cell] += rDPtr[cell]*bc[face]*recvBuf[face];
            }
        }
        Pout<< "** after higherNbrs_:" << flatOutput(wA) << endl;

        // Do non-processor boundaries
        for (label colouri = 0; colouri < nColours_; colouri++)
        {
            Pout<< "** starting backwards colour:" << colouri << endl;

            check();
            for (const label inti : higherGlobal_)
            {
                // Note: interfaceBouCoeffs is -upperPtr
                const auto& intf = interfaces[inti].interface();
                auto& coeff = coeffs_[inti];
                const auto& faceCells = intf.faceCells();
                const auto& bc = interfaceBouCoeffs[inti];
                forAll(coeff, face)
                {
                    const label cell = faceCells[face];
                    if (cellColour_[cell] == colouri)
                    {
                        coeff[face] = rDPtr[cell]*bc[face];
                    }
                }
            }

            updateMatrixInterfaces
            (
                false,           // add to rD
                coeffs_,
                higherGlobal_,  //global_,
                wA,             // send data
                wA,             // result
                0               // cmpt
            );
            // Reset value to zero
            zeroCoeffs(higherGlobal_);
            Pout<< "** after higherGlobal_:" << flatOutput(wA) << endl;

            for (label face=nFacesM1; face>=0; face--)
            {
                const label cell = uPtr[face];
                if (cellColour_[cell] == colouri)
                {
                    wAPtr[lPtr[face]] -=
                        rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
                }
            }

            Pout<< "** finished backwards colour:" << colouri << nl << endl;
        }
    }


    Pout<< "** after backwards sweep:" << flatOutput(wA) << endl;

    // Start writes of wA (using sendBufs)
    send(lowerNbrs_, wA);
    UPstream::waitRequests(startOfRequests);
}


// ************************************************************************* //
