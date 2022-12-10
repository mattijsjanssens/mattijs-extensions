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

void Foam::distributedDICPreconditioner2::walkFront
(
    const lduMesh& lm,
    DynamicList<label>& front,
    labelList& cellColour
) const
{
    // Colour with the (min) coupled (global) patch

    const lduAddressing& addr = lm.lduAddr();
    const label* const __restrict__ uPtr = addr.upperAddr().begin();
    const label* const __restrict__ lPtr = addr.lowerAddr().begin();
    const label* const __restrict__ ownStartPtr = addr.ownerStartAddr().begin();
    const label* const __restrict__ losortStartAddrPtr =
        addr.losortStartAddr().begin();
    const label* const __restrict__ losortAddrPtr = addr.losortAddr().begin();

    DynamicList<label> newFront;
    while (true)
    {
        newFront.clear();
        for (const label celli : front)
        {
            const label colouri = cellColour[celli];

            {
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
                        cellColour[nbr] = colouri;
                        newFront.append(nbr);
                    }
                }
            }
            {
                const label fStart = losortStartAddrPtr[celli];
                const label fEnd = losortStartAddrPtr[celli + 1];

                for (label i=fStart; i<fEnd; i++)
                {
                    label facei = losortAddrPtr[i];
                    const label nbr =
                    (
                        lPtr[facei] == celli
                      ? uPtr[facei]
                      : lPtr[facei]
                    );
                    if (cellColour[nbr] == -1)
                    {
                        cellColour[nbr] = colouri;
                        newFront.append(nbr);
                    }
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


Foam::label Foam::distributedDICPreconditioner2::meshColour
(
    const lduMesh& lm,
    labelList& cellColour,
    labelList& patchToColour
) const
{
    // Colour with the (min) coupled (global) patch. Compact to have colour 0
    // if a single region. Problematic if same patch on multiple disconnected
    // regions.

    const lduAddressing& addr = lm.lduAddr();
    const auto& interfaces = solver_.interfaces();

    const label nCells = addr.size();

    patchToColour.setSize(interfaces.size());
    patchToColour = -1;
    cellColour.setSize(nCells);
    cellColour = -1;


    label colouri = 0;

    // Starting front from patch faceCells
    DynamicList<label> front;
    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            const auto& intf = interfaces[inti].interface();
            if (!isA<const processorLduInterface>(intf))
            {
                // 'global' interface. Seed faceCells with patch index

                if (patchToColour[inti] == -1)
                {
                    patchToColour[inti] = colouri++;
                }

                const auto& fc = intf.faceCells();
                forAll(fc, face)
                {
                    const label cell = fc[face];
                    if (cellColour[cell] != patchToColour[inti])
                    {
                        cellColour[cell] = patchToColour[inti];
                        front.append(cell);
                    }
                }
            }
        }
    }


    // Walk current mesh
    walkFront(lm, front, cellColour);


    // Any unvisited cell is still labelMax. Put into colour 0.
    for (auto& colour : cellColour)
    {
        if (colour == -1)
        {
            colour = 0;
        }
    }
    return colouri;
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
    UPstream::msgType() += 321;

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

    UPstream::msgType() -= 321;
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
            ppp->tag()+70,          // random offset to avoid comms conflict
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

        Pout<< "Dest proc:" << ppp->neighbProcNo()
            << " interface:" << inti
            << " buf:" << flatOutput(sendBuf)
            << endl;


        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            ppp->neighbProcNo(),
            sendBuf.cdata_bytes(),
            sendBuf.size_bytes(),
            ppp->tag()+70,
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

    if (selectedInterfaces.size())
    {
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
            selectedInterfaces,
            psi,                // send data
            psi,                // result
            0                   // cmpt
        );

        auto& colourBuf = colourBufs_[colouri];
        colourBuf.setSize(interfaces.size());
        for (const label inti : selectedInterfaces)
        {
            const auto& intf = interfaces[inti].interface();
            const auto& fc = intf.faceCells();
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
                << " buf:" << flatOutput(recvBuf)
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


    // - divide the internal mesh into domains/colour, similar to domain
    //   decomposition
    // - we could use subsetted bits of the mesh or just loop over only
    //   the cells of the current domain
    // - a domain only uses boundary values of lower numbered domains
    // - this also limits the interfaces that need to be evaluated since
    //   we assume an interface only changes local faceCells and these
    //   are all of the same colour

    for (label colouri = 0; colouri < nColours_; colouri++)
    {
        Pout<< "** starting colour:" << colouri << endl;

        if (cellColourPtr_.valid() && lowerGlobalRecv_[colouri].size())
        {
            const auto& cellColour = cellColourPtr_();
            for (const label inti : lowerGlobalRecv_[colouri])
            {
                const auto& intf = interfaces[inti].interface();
                // TBD: do not use patch faceCells but passed-in addressing?
                const auto& faceCells = intf.faceCells();
                const auto& recvBuf = colourBufs_[colouri][inti];
                const auto& bc = interfaceBouCoeffs[inti];

                Pout<< "    for interface:" << inti
                    << " retrieving slot:" << colouri
                    << " buf:" << recvBuf
                    << " inserting into cells:" << flatOutput(faceCells)
                    << " with coeffs:" << flatOutput((bc*bc)())
                    << endl;

                forAll(recvBuf, face)
                {
                    // Note:interfaceBouCoeffs is -upperPtr
                    const label cell = faceCells[face];
                    if (cellColour[cell] != colouri)
                    {
                        FatalErrorInFunction << " problem" << exit(FatalError);
                    }
                    rD[cell] -= bc[face]*bc[face]/recvBuf[face];
                }
            }
            Pout<< "** after lowerGlobalRecv_:" << flatOutput(rD) << endl;
        }


        const label nFaces = matrix.upper().size();
        if (cellColourPtr_.valid())
        {
            const auto& cellColour = cellColourPtr_();
            for (label face=0; face<nFaces; face++)
            {
                const label cell = lPtr[face];
                if (cellColour[cell] == colouri)
                {
                    rD[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rD[cell];
                }
            }
        }
        else
        {
            for (label face=0; face<nFaces; face++)
            {
                rD[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rD[lPtr[face]];
            }
        }

        // Store effect of exchanging rD to higher interfaces in colourBufs_
        if (cellColourPtr_.valid())
        {
            sendGlobal(higherGlobalSend_[colouri], rD, higherColour_[colouri]);
        }

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
    const dictionary& dict
)
:
    lduMatrix::preconditioner(sol),
    coupled_(dict.getOrDefault<bool>("coupled", true, keyType::LITERAL)),
    rD_(sol.matrix().diag().size())
{
    const lduMesh& mesh = sol.matrix().mesh();
    const auto& interfaces = sol.interfaces();
    const auto& interfaceBouCoeffs = sol.interfaceBouCoeffs();

    // Allocate buffers
    // ~~~~~~~~~~~~~~~~

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
    //const label myColour = colours[Pstream::myProcNo()];
    const label myColour = Pstream::myProcNo();

    bool haveCyclicAMI = false;
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
                if (nbrColour < myColour)
                {
                    lowerNbrs_.append(inti);
                }
                else if (nbrColour > myColour)
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
            else if (isA<const cyclicAMILduInterface>(intf))
            {
                haveCyclicAMI = true;
            }   
        }
    }


    // Determine local colouring/zoneing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Start off with all cells colour 0
    nColours_ = 1;
    cellColourPtr_.clear();
    if (coupled_ && haveCyclicAMI)
    {
        labelList patchToColour;
        cellColourPtr_.reset(new labelList(0));
        nColours_ = meshColour(mesh, cellColourPtr_(), patchToColour);

        Pout<< "Determining nColours:" << nColours_
            << " cellColour:" << flatOutput(cellColourPtr_()) << endl;

        lowerGlobalRecv_.setSize(nColours_);
        lowerGlobalSend_.setSize(nColours_);
        lowerColour_.setSize(nColours_, -1);
        higherGlobalRecv_.setSize(nColours_);
        higherGlobalSend_.setSize(nColours_);
        higherColour_.setSize(nColours_, -1);
        colourBufs_.setSize(nColours_);

        forAll(interfaces, inti)
        {
            if (interfaces.set(inti))
            {
                const auto& intf = interfaces[inti].interface();
                const auto* AMIpp = isA<const cyclicAMILduInterface>(intf);
                if (AMIpp)
                {
                    const label colouri = patchToColour[inti];
                    const label nbrInti = AMIpp->neighbPatchID();
                    const label nbrColouri = patchToColour[nbrInti];
                    const auto& bc = interfaceBouCoeffs[inti];

                    if (!colourBufs_.set(colouri))
                    {
                        colourBufs_.set
                        (
                            colouri,
                            new FieldField<Field, solveScalar>
                            (
                                interfaces.size()
                            )
                        );
                    }
                    auto& colourBuf = colourBufs_[colouri];
                    if (!colourBuf.set(inti))
                    {
                        colourBuf.set
                        (
                            inti,
                            new solveScalarField(bc.size(), Zero)
                        );
                    }

                    if (colouri < nbrColouri)
                    {
                        // Send to higher
                        higherGlobalSend_[colouri].append(nbrInti);
                        higherColour_[colouri] = nbrColouri;
                        // Receive from higher
                        higherGlobalRecv_[colouri].append(inti);
                    }
                    else
                    {
                        // Send to lower
                        lowerGlobalSend_[colouri].append(nbrInti);
                        lowerColour_[colouri] = nbrColouri;
                        // Receive from lower
                        lowerGlobalRecv_[colouri].append(inti);
                    }
                }
            }
        }

        DebugVar(lowerGlobalRecv_);
        DebugVar(lowerColour_);
        DebugVar(higherGlobalSend_);
        DebugVar(higherColour_);
    }

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

            if (cellColourPtr_.valid() && lowerGlobalRecv_[colouri].size())
            {
                const auto& cellColour = cellColourPtr_();
                for (const label inti : lowerGlobalRecv_[colouri])
                {
                    Pout<< "    for interface:" << inti
                        << " retrieving slot:" << colouri << endl;

                    const auto& intf = interfaces[inti].interface();
                    const auto& faceCells = intf.faceCells();
                    const auto& recvBuf = colourBufs_[colouri][inti];
                    const auto& bc = interfaceBouCoeffs[inti];

                    Pout<< "    for interface:" << inti
                        << " retrieving slot:" << colouri
                        << " buf:" << recvBuf
                        << " inserting into cells:" << flatOutput(faceCells)
                        << " with coeffs:" << flatOutput(bc)
                        << endl;

                    forAll(recvBuf, face)
                    {
                        // Note:interfaceBouCoeffs is -upperPtr
                        const label cell = faceCells[face];
                        if (cellColour[cell] != colouri)
                        {
                            FatalErrorInFunction
                                << " problem" << exit(FatalError);
                        }
                        wAPtr[cell] += rDPtr[cell]*bc[face]*recvBuf[face];
                    }
                }
                Pout<< "** after lowerGlobalRecv_:" << flatOutput(wA) << endl;
            }

            if (cellColourPtr_.valid())
            {
                const auto& cellColour = cellColourPtr_();
                for (label face=0; face<nFaces; face++)
                {
                    const label cell = lPtr[face];
                    if (cellColour[cell] == colouri)
                    {
                        wAPtr[uPtr[face]] -=
                            rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[cell];
                    }
                }
            }
            else
            {
                for (label face=0; face<nFaces; face++)
                {
                    wAPtr[uPtr[face]] -=
                        rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
                }
            }

            // Store effect of exchanging rD to higher interfaces
            // in colourBufs_
            if (cellColourPtr_.valid())
            {
                sendGlobal
                (
                    higherGlobalSend_[colouri],
                    wA,
                    higherColour_[colouri]
                );
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

            Pout<< "inti:" << inti
                << " buf:" << flatOutput(recvBuf)
                << " inserting into cells:" << flatOutput(faceCells)
                << " with coeffs:" << flatOutput(bc)
                << endl;

            forAll(recvBuf, face)
            {
                // Note: interfaceBouCoeffs is -upperPtr
                const label cell = faceCells[face];
                wAPtr[cell] += rDPtr[cell]*bc[face]*recvBuf[face];
            }
        }
        Pout<< "** after higherNbrs_:" << flatOutput(wA) << endl;

        // Do non-processor boundaries
        for (label colouri = nColours_-1; colouri >= 0; colouri--)
        {
            Pout<< "** starting backwards colour:" << colouri << endl;

            if (cellColourPtr_.valid() && higherGlobalRecv_[colouri].size())
            {
                const auto& cellColour = cellColourPtr_();
                for (const label inti : higherGlobalRecv_[colouri])
                {
                    const auto& intf = interfaces[inti].interface();
                    const auto& faceCells = intf.faceCells();
                    const auto& recvBuf = colourBufs_[colouri][inti];
                    const auto& bc = interfaceBouCoeffs[inti];

                    Pout<< "    for interface:" << inti
                        << " retrieving slot:" << colouri
                        << " buf:" << recvBuf
                        << " inserting into cells:" << flatOutput(faceCells)
                        << " with coeffs:" << flatOutput(bc)
                        << endl;
                    forAll(recvBuf, face)
                    {
                        // Note:interfaceBouCoeffs is -upperPtr
                        const label cell = faceCells[face];
                        if (cellColour[cell] != colouri)
                        {
                            FatalErrorInFunction
                                << " problem" << exit(FatalError);
                        }
                        wAPtr[cell] += rDPtr[cell]*bc[face]*recvBuf[face];
                    }
                }
                Pout<< "** after higherGlobalRecv_:" << flatOutput(wA) << endl;
            }

            if (cellColourPtr_.valid())
            {
                const auto& cellColour = cellColourPtr_();
                for (label face=nFacesM1; face>=0; face--)
                {
                    const label cell = uPtr[face];
                    if (cellColour[cell] == colouri)
                    {
                        // Note: lower cell guaranteed in same colour
                        wAPtr[lPtr[face]] -=
                            rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[cell];
                    }
                }
            }
            else
            {
                for (label face=nFacesM1; face>=0; face--)
                {
                    wAPtr[lPtr[face]] -=
                        rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
                }
            }

            // Store effect of exchanging rD to higher interfaces
            // in colourBufs_
            if (cellColourPtr_.valid())
            {
                sendGlobal
                (
                    lowerGlobalSend_[colouri],
                    wA,
                    lowerColour_[colouri]
                );
            }

            Pout<< "** finished backwards colour:" << colouri << nl << endl;
        }
    }

    Pout<< "** after backwards sweep:" << flatOutput(wA) << endl;

    // Start writes of wA (using sendBufs)
    send(lowerNbrs_, wA);
    //UPstream::waitRequests(startOfRequests);

    Pout<< nl << endl;

}


// ************************************************************************* //
