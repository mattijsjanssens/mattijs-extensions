/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "distributedDICPreconditioner-Request.H"
#include "processorLduInterface.H"
#include "processorColour.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedDICPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<distributedDICPreconditioner>
        adddistributedDICPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distributedDICPreconditioner::receive
(
    const labelList& selectedInterfaces,
    DynamicList<UPstream::Request>& requests
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
        recvBuf.resize_nocopy(interfaceBouCoeffs[inti].size());

        requests.push_back(UPstream::Request());
        UIPstream::read
        (
            requests.back(),
            ppp->neighbProcNo(),
            recvBuf.data_bytes(),
            recvBuf.size_bytes(),
            ppp->tag()+70,          // random offset
            ppp->comm()
        );
    }
}


void Foam::distributedDICPreconditioner::send
(
    const labelList& selectedInterfaces,
    const solveScalarField& psiInternal,
    DynamicList<UPstream::Request>& requests
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

        sendBuf.resize_nocopy(faceCells.size());
        forAll(faceCells, face)
        {
            sendBuf[face] = psiInternal[faceCells[face]];
        }

        requests.push_back(UPstream::Request());
        UOPstream::write
        (
            requests.back(),
            ppp->neighbProcNo(),
            sendBuf.cdata_bytes(),
            sendBuf.size_bytes(),
            ppp->tag()+70,          // random offset
            ppp->comm()
        );
    }
}


void Foam::distributedDICPreconditioner::wait
(
    DynamicList<UPstream::Request>& requests
) const
{
    UPstream::waitRequests(requests);
    requests.clear();
}


void Foam::distributedDICPreconditioner::calcReciprocalD
(
    solveScalarField& rD
) const
{
    const auto& interfaces = solver_.interfaces();
    const auto& interfaceBouCoeffs = solver_.interfaceBouCoeffs();
    const auto& matrix = solver_.matrix();
    const auto& lduAddr = matrix.lduAddr();

    // Make sure no outstanding receives
    wait(lowerRecvRequests_);

    // Start reads (into recvBufs_)
    receive(lowerNbrs_, lowerRecvRequests_);

    // Start with diagonal
    const scalarField& diag = matrix.diag();
    rD.resize_nocopy(diag.size());
    std::copy(diag.begin(), diag.end(), rD.begin());


    const label* const __restrict__ uPtr = lduAddr.upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr.lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    // Subtract coupled contributions
    {
        // Wait for finish. Received result in recvBufs
        wait(lowerRecvRequests_);

        for (const label inti : lowerNbrs_)
        {
            const auto& intf = interfaces[inti].interface();
            // TBD: do not use patch faceCells but passed-in addressing?
            const auto& faceCells = intf.faceCells();
            const auto& recvBuf = recvBufs_[inti];
            const auto& bc = interfaceBouCoeffs[inti];

            forAll(recvBuf, face)
            {
                // Note:interfaceBouCoeffs is -upperPtr
                rD[faceCells[face]] -= bc[face]*bc[face]/recvBuf[face];
            }
        }
    }


    const label nFaces = matrix.upper().size();
    for (label face=0; face<nFaces; face++)
    {
        rD[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rD[lPtr[face]];
    }


    // Make sure no outstanding sends
    wait(higherSendRequests_);


    // Start writes of rD (using sendBufs)
    send(higherNbrs_, rD, higherSendRequests_);


    // Calculate the reciprocal of the preconditioned diagonal
    const label nCells = rD.size();

    for (label cell=0; cell<nCells; cell++)
    {
        rD[cell] = 1.0/rD[cell];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedDICPreconditioner::distributedDICPreconditioner
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

    sendBufs_.setSize(interfaces.size());
    recvBufs_.setSize(interfaces.size());
    lowerNbrs_.setCapacity(interfaces.size());
    lowerSendRequests_.setCapacity(interfaces.size());
    lowerRecvRequests_.setCapacity(interfaces.size());
    higherNbrs_.setCapacity(interfaces.size());
    higherSendRequests_.setCapacity(interfaces.size());
    higherRecvRequests_.setCapacity(interfaces.size());
    forAll(interfaceBouCoeffs, inti)
    {
        if (interfaces.set(inti))
        {
            const auto& intf = interfaces[inti].interface();
            const auto* ppp = isA<const processorLduInterface>(intf);
            if (ppp)
            {
                const auto& bc = interfaceBouCoeffs[inti];
                sendBufs_.set(inti, new scalarField(bc.size(), Zero));
                recvBufs_.set(inti, new scalarField(bc.size(), Zero));
            }
        }
    }


    // Determine processor colouring
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //const processorColour& colours = processorColour::New(mesh);
    //const label nColours = 
    processorColour::colour(mesh, colours_);
    //Pout<< "Calculated colouring with " << nColours << endl;
    const label colouri = colours_[Pstream::myProcNo(mesh.comm())];
    //const label colouri = Pstream::myProcNo(mesh.comm());

    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            const auto& intf = interfaces[inti].interface();
            const auto* ppp = isA<const processorLduInterface>(intf);
            if (ppp)
            {
                const label nbrColour = colours_[ppp->neighbProcNo()];
                //const label nbrColour = ppp->neighbProcNo();
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
        }
    }

    calcReciprocalD(rD_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributedDICPreconditioner::~distributedDICPreconditioner()
{
    // Make sure no requests are still in flight (since preconditioner
    // ends with send-to-lower)

    wait(lowerSendRequests_);
    wait(higherSendRequests_);
    wait(lowerRecvRequests_);
    wait(higherRecvRequests_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distributedDICPreconditioner::precondition
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

    // Make sure no receives are still in flight (should not happen)
    wait(lowerRecvRequests_);

    // Start reads (into recvBufs)
    receive(lowerNbrs_, lowerRecvRequests_);

    // Initialise 'internal' cells
    for (label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }

    // Do 'halo' contributions from lower numbered procs
    {
        // Wait for finish. Received result in recvBufs
        wait(lowerRecvRequests_);

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
    }

    for (label face=0; face<nFaces; face++)
    {
        wAPtr[uPtr[face]] -= rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
    }

    // Make sure no outstanding sends from previous iteration
    wait(higherSendRequests_);

    // Start writes of wA (using sendBufs)
    send(higherNbrs_, wA, higherSendRequests_);


    // Backward sweep
    // ~~~~~~~~~~~~~~

    // Make sure no outstanding receives
    wait(higherRecvRequests_);

    // Start receives
    receive(higherNbrs_, higherRecvRequests_);

    {
        // Wait until receives finished
        wait(higherRecvRequests_);

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
    }


    for (label face=nFacesM1; face>=0; face--)
    {
        wAPtr[lPtr[face]] -= rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
    }

    // Make sure no outstanding sends
    wait(lowerSendRequests_);

    // Start writes of wA (using sendBufs)
    send(lowerNbrs_, wA, lowerSendRequests_);
}


// ************************************************************************* //
