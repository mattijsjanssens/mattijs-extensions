/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2022 M. Janssens
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

#include "DICPCG.H"
#include "processorLduInterface.H"
#include "processorColour.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DICPCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<DICPCG>
        addDICPCGSymMatrixConstructorToTable_;
}


void Foam::DICPCG::receive
(
    const labelList& selectedInterfaces,
    FieldField<Field, scalar>& recvBufs,    // Receive buffer
    label& startOfRequests
) const
{
    if (startOfRequests == -1)
    {
        startOfRequests = Pstream::nRequests();
    }

    // Start reads
    for (const label inti : selectedInterfaces)
    {
        const auto& intf = interfaces_[inti].interface();
        const auto* ppp = isA<const processorLduInterface>(intf);

        scalarField& recvBuf = recvBufs[inti];
        recvBuf.setSize(interfaceBouCoeffs_[inti].size());

        Pout<< "** starting read at interface:" << inti
            << " from " << ppp->neighbProcNo() << endl;

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


void Foam::DICPCG::send
(
    const labelList& selectedInterfaces,
    const solveScalarField& psiInternal,
    FieldField<Field, scalar>& sendBufs,
    label& startOfRequests
) const
{
    if (startOfRequests == -1)
    {
        startOfRequests = Pstream::nRequests();
    }

    // Start writes
    for (const label inti : selectedInterfaces)
    {
        const auto& intf = interfaces_[inti].interface();
        const auto* ppp = isA<const processorLduInterface>(intf);
        const auto& faceCells = intf.faceCells();

        scalarField& sendBuf = sendBufs[inti];

        Pout<< "** starting send at interface:" << inti
            << " to " << ppp->neighbProcNo() << endl;

        sendBuf.setSize(faceCells.size());
        forAll(faceCells, face)
        {
            sendBuf[face] = psiInternal[faceCells[face]];
        }

        Pout<< "at int:" << inti
            << " sending to proc:" << ppp->neighbProcNo()
            << " actual data:" << flatOutput(sendBuf)
            << endl;

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


void Foam::DICPCG::calcReciprocalD
(
    solveScalarField& rD,
    const lduMatrix& matrix,
    const labelList& lowerInterfaces,
    const labelList& higherInterfaces,
    FieldField<Field, scalar>& sendBufs,    // Send buffer
    FieldField<Field, scalar>& recvBufs,    // Receive buffer
    label& startOfRequests
) const
{
    // Start swapping remote contributions
    startOfRequests = Pstream::nRequests();

    // Start reads (into recvBufs)
    Pout<< "Starting read from interfaces:" << flatOutput(lowerInterfaces)
        << endl;
    receive(lowerInterfaces, recvBufs, startOfRequests);


    const scalarField& diag = matrix.diag();
    rD.setSize(diag.size());
    std::copy(diag.begin(), diag.end(), rD.begin());


    const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    // Subtract coupled contributions
    if (startOfRequests != -1)
    {
        // Wait for finish. Received result in recvBufs
        UPstream::waitRequests(startOfRequests);
        startOfRequests = -1;
    }

    for (const label inti : lowerInterfaces)
    {
        const auto& intf = interfaces_[inti].interface();
        // TBD: do not use patch faceCells but passed-in addressing?
        const auto& faceCells = intf.faceCells();
        const auto& recvBuf = recvBufs[inti];
        const auto& bc = interfaceBouCoeffs_[inti];

        const auto* ppp = isA<const processorLduInterface>(intf);
        Pout<< "at int:" << inti
            << " received from proc:" << ppp->neighbProcNo()
            << " actual data:" << flatOutput(recvBuf)
            << endl;

        forAll(recvBuf, face)
        {
            // Note:interfaceBouCoeffs_ is -upperPtr
            rD[faceCells[face]] -= bc[face]*bc[face]/recvBuf[face];
        }
    }


    const label nFaces = matrix.upper().size();
    for (label face=0; face<nFaces; face++)
    {
        rD[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rD[lPtr[face]];
    }


    // Start writes of rD (using sendBufs)
    send(higherInterfaces, rD, sendBufs, startOfRequests);

    // Calculate the reciprocal of the preconditioned diagonal
    const label nCells = rD.size();

    for (label cell=0; cell<nCells; cell++)
    {
        rD[cell] = 1.0/rD[cell];
    }

    if (startOfRequests != -1)
    {
        // Wait for finish
        UPstream::waitRequests(startOfRequests);
        startOfRequests = -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DICPCG::DICPCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::DICPCG::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    label nCells = psi.size();

    scalar* __restrict__ psiPtr = psi.begin();

    scalarField pA(nCells);
    scalar* __restrict__ pAPtr = pA.begin();

    scalarField wA(nCells);
    scalar* __restrict__ wAPtr = wA.begin();

    scalar wArA = solverPerf.great_;
    scalar wArAold = wArA;

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    scalarField rA(source - wA);
    scalar* __restrict__ rAPtr = rA.begin();

    // --- Calculate normalisation factor
    scalar normFactor = this->normFactor(psi, source, wA, pA);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();


    // --- Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
        //// --- Select and construct the preconditioner
        //autoPtr<lduMatrix::preconditioner> preconPtr =
        //lduMatrix::preconditioner::New
        //(
        //    *this,
        //    controlDict_
        //);

//XXXXX
        // Work array
        FieldField<Field, scalar> sendBufs(interfaceBouCoeffs_.size());
        FieldField<Field, scalar> recvBufs(interfaceBouCoeffs_.size());
        forAll(interfaceBouCoeffs_, inti)
        {
            if (interfaces_.set(inti))
            {
                const auto& bc = interfaceBouCoeffs_[inti];
                sendBufs.set(inti, new scalarField(bc.size(), Zero));
                recvBufs.set(inti, new scalarField(bc.size(), Zero));
            }
        }
        label startOfRequests = -1;

        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
            //preconPtr->precondition(wA, rA, cmpt);
//XXXXXX
            {
                // TBD: replace with proper colouring
                //const labelList colours(identity(Pstream::nProcs()));
                //const label nColours = max(colours)+1;
                const lduMesh& mesh = matrix().mesh();
                const processorColour& colours = processorColour::New(mesh);
                //const label nColours = colours.nColours();
                const label colouri = colours[Pstream::myProcNo()];

                // Select interfaces involving colouri
                // - me to lower-coloured procs (= already solved for)
                // - lower-coloured procs to me
                DynamicList<label> lowerNbrs;
                DynamicList<label> higherNbrs;
                forAll(interfaceBouCoeffs_, inti)
                {
                    if (interfaces_.set(inti))
                    {
                        const auto& intf = interfaces_[inti].interface();
                        const auto* ppp =
                            isA<const processorLduInterface>(intf);
                        if (ppp)
                        {
                            if (colours[ppp->neighbProcNo()] < colouri)
                            {
                                // Neighbour not yet updated
                                lowerNbrs.append(inti);
                            }
                            else if (colours[ppp->neighbProcNo()] > colouri)
                            {
                                higherNbrs.append(inti);
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

                // Include not-solved for neighbours
                Pout<< "** Diagonal" << endl;
                scalarField rD_;
                //calcReciprocalD(rD_, matrix(), lowerNbrs, coeffs, cmpt);
                {
                    calcReciprocalD
                    (
                        rD_,
                        matrix(),
                        lowerNbrs,
                        higherNbrs,
                        sendBufs,    // Send buffer
                        recvBufs,    // Receive buffer
                        startOfRequests
                    );
                    DebugVar(rD_);
                }



                scalar* __restrict__ wAPtr = wA.begin();
                const scalar* __restrict__ rAPtr = rA.begin();
                const scalar* __restrict__ rDPtr = rD_.begin();

                const label* const __restrict__ uPtr =
                    matrix().lduAddr().upperAddr().begin();
                const label* const __restrict__ lPtr =
                    matrix().lduAddr().lowerAddr().begin();
                const scalar* const __restrict__ upperPtr =
                    matrix().upper().begin();

                label nCells = wA.size();
                label nFaces = matrix().upper().size();
                label nFacesM1 = nFaces - 1;


                label startOfRequests = Pstream::nRequests();

                // Start reads (into recvBufs)
                Pout<< "Starting read from interfaces:"
                    << flatOutput(lowerNbrs) << endl;
                receive(lowerNbrs, recvBufs, startOfRequests);
                UPstream::waitRequests(startOfRequests);
                startOfRequests = -1;

                Pout<< "** forwards sweep" << endl;

                // Initialise 'internal' cells
                for (label cell=0; cell<nCells; cell++)
                {
                    wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
                }

                // Do 'halo' contributions from lower numbered procs
                if (startOfRequests != -1)
                {
                    // Wait for finish. Received result in recvBufs
                    UPstream::waitRequests(startOfRequests);
                    startOfRequests = -1;
                }

                for (const label inti : lowerNbrs)
                {
                    const auto& intf = interfaces_[inti].interface();
                    // TBD: do not use patch faceCells but passed-in
                    // addressing?
                    const auto& faceCells = intf.faceCells();
                    const auto& recvBuf = recvBufs[inti];
                    const auto& bc = interfaceBouCoeffs_[inti];

                    const auto* ppp =
                        isA<const processorLduInterface>(intf);
                    Pout<< "at int:" << inti
                        << " received from proc:" << ppp->neighbProcNo()
                        << " actual data:" << flatOutput(recvBuf)
                        << endl;

                    forAll(recvBuf, face)
                    {
                        // Note: interfaceBouCoeffs_ is -upperPtr
                        const label cell = faceCells[face];

                        Pout<< "From remote:" << rDPtr[cell]*bc[face]*recvBuf[face]
                            << endl;

                        wAPtr[cell] += rDPtr[cell]*bc[face]*recvBuf[face];
                    }
                }

                // Do 'internal' faces, forward sweep
                for (label face=0; face<nFaces; face++)
                {
                    Pout<< "From lower cell:"
                        << -rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]]
                        << endl;

                    wAPtr[uPtr[face]] -=
                        rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
                }

                // Start send/receives
                {
                    // Start writes of wA (using sendBufs)
                    send(higherNbrs, wA, sendBufs, startOfRequests);
                    UPstream::waitRequests(startOfRequests);
                    startOfRequests = -1;
                }


                Pout<< "** backwards sweep" << endl;

                // Do 'halo' contributions from higher numbered procs
                if (startOfRequests != -1)
                {
                    // Wait for finish. Received result in recvBufs
                    UPstream::waitRequests(startOfRequests);
                    startOfRequests = -1;
                }
                for (const label inti : higherNbrs)
                {
                    const auto& intf = interfaces_[inti].interface();
                    // TBD: do not use patch faceCells but passed-in
                    // addressing?
                    const auto& faceCells = intf.faceCells();
                    const auto& recvBuf = recvBufs[inti];
                    const auto& bc = interfaceBouCoeffs_[inti];

                    const auto* ppp =
                        isA<const processorLduInterface>(intf);
                    Pout<< "at int:" << inti
                        << " received from proc:" << ppp->neighbProcNo()
                        << " actual data:" << flatOutput(recvBuf)
                        << endl;

                    forAll(recvBuf, face)
                    {
                        // Note: interfaceBouCoeffs_ is -upperPtr
                        const label cell = faceCells[face];

                        Pout<< "From remote:" << rDPtr[cell]*bc[face]*recvBuf[face]
                            << endl;

                        wAPtr[cell] += rDPtr[cell]*bc[face]*recvBuf[face];
                    }
                }

                // Do 'internal' faces, backwards sweep
                for (label face=nFacesM1; face>=0; face--)
                {
                    Pout<< "From upper cell:"
                        << -rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]]
                        << endl;

                    wAPtr[lPtr[face]] -=
                        rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
                }

                {
                    // Start writes of wA (using sendBufs)
                    send(lowerNbrs, wA, sendBufs, startOfRequests);
                    UPstream::waitRequests(startOfRequests);
                    startOfRequests = -1;
                }
                Pout<< endl;
            }
//XXXXXX

            // --- Update search directions:
            wArA = gSumProd(wA, rA, matrix().mesh().comm());

            if (solverPerf.nIterations() == 0)
            {
                for (label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell];
                }
            }
            else
            {
                scalar beta = wArA/wArAold;

                for (label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                }
            }


            // --- Update preconditioned residual
            matrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);

            scalar wApA = gSumProd(wA, pA, matrix().mesh().comm());


            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;


            // --- Update solution and residual:

            scalar alpha = wArA/wApA;

            for (label cell=0; cell<nCells; cell++)
            {
                psiPtr[cell] += alpha*pAPtr[cell];
                rAPtr[cell] -= alpha*wAPtr[cell];
            }

            solverPerf.finalResidual() =
                gSumMag(rA, matrix().mesh().comm())
               /normFactor;

        } while
        (
            (
                solverPerf.nIterations()++ < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    return solverPerf;
}


// ************************************************************************* //
