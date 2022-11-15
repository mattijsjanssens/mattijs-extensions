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

    See chapter 12.6 / page 374 of Iterative Methods for Sparse Linear Systems


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


void Foam::DICPCG::calcReciprocalD
(
    solveScalarField& rD,
    const lduMatrix& matrix
) const
{
    solveScalar* __restrict__ rDPtr = rD.begin();

    const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    // Calculate the DIC diagonal
    const label nFaces = matrix.upper().size();
    for (label face=0; face<nFaces; face++)
    {
        const scalar s = upperPtr[face]*upperPtr[face]/rDPtr[lPtr[face]];

//        Pout<< "calcReciprocalD : For face:" << face
//            << " adapting diagonal for uppercell:" << uPtr[face]
//            << " weight:" << upperPtr[face]
//            << " contributions from lowercell:" << lPtr[face]
//            << " lowerContrib:" << s
//            << endl;

        rDPtr[uPtr[face]] -= s;
    }


    // Calculate the reciprocal of the preconditioned diagonal
    const label nCells = rD.size();

    for (label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }
}
void Foam::DICPCG::calcReciprocalD
(
    solveScalarField& rD,
    const lduMatrix& matrix,
    const labelList& selectedInterfaces,
    FieldField<Field, scalar>& coeffs,  // work
    const direction cmpt
) const
{
    const scalarField& diag = matrix.diag();
    rD.setSize(diag.size());
    std::copy(diag.begin(), diag.end(), rD.begin());


    // Subtract coupled contributions
    {
        scalarField invRd(matrix.diag());
        forAll(invRd, cell)
        {
            invRd[cell] = 1.0/invRd[cell];
        }


        forAll(interfaces_, inti)
        {
            if (interfaces_.set(inti))
            {
                coeffs[inti] = Zero;
            }
        }
        for (const label inti : selectedInterfaces)
        {
            const auto& bc = interfaceBouCoeffs_[inti];
            scalarField& coeff = coeffs[inti];

            forAll(bc, i)
            {
                coeff[i] = bc[i]*bc[i];
            }
        }

        const label startRequest = Pstream::nRequests();
        matrix.initMatrixInterfaces
        (
            true,   // subtract remote contributions
            coeffs,
            interfaces_,
            invRd,
            rD,
            cmpt                
        );
        matrix.updateMatrixInterfaces
        (
            true,   // subtract remote contributions
            coeffs,
            interfaces_,
            invRd,
            rD,
            cmpt,
            startRequest              
        );
    }

    const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    {
        const label nFaces = matrix.upper().size();

        for (label face=0; face<nFaces; face++)
        {
            rD[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rD[lPtr[face]];
        }


        // Calculate the reciprocal of the preconditioned diagonal
        const label nCells = rD.size();

        for (label cell=0; cell<nCells; cell++)
        {
            rD[cell] = 1.0/rD[cell];
        }
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
{
//    const scalarField& diag = matrix.diag();
//
//    scalarField rD(diag.size());
//    std::copy(diag.begin(), diag.end(), rD.begin());
//
//    calcReciprocalD(rD, matrix);
//
//    Pout<< "** non-parallel rD:" << flatOutput(rD) << endl;
}


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
        //{
        //    scalarField nonParRD(matrix().diag());
        //    calcReciprocalD(nonParRD, matrix());
        //    DebugVar(nonParRD);
        //}
//XXXX



        // Work array
        FieldField<Field, scalar> coeffs(interfaceBouCoeffs_.size());
        forAll(interfaceBouCoeffs_, inti)
        {
            if (interfaces_.set(inti))
            {
                const auto& bc = interfaceBouCoeffs_[inti];
                coeffs.set(inti, new scalarField(bc.size(), Zero));
            }
        }



        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
            //preconPtr->precondition(wA, rA, cmpt);

//XXXXXX
            // TBD: replace with proper colouring
            const lduMesh& mesh = matrix().mesh();
            const processorColour& colours = processorColour::New(mesh);
            const label nColours = colours.nColours();

            for (label colouri = 0; colouri < nColours; colouri++)
            {
                // Select interfaces involving colouri
                // - me to lower-coloured procs (= already solved for)
                // - lower-coloured procs to me
                DynamicList<label> meFromHigherNbrs;
                DynamicList<label> meFromLowerNbrs;
                DynamicList<label> higherFromMeNbrs;
                forAll(interfaceBouCoeffs_, inti)
                {
                    if (interfaces_.set(inti))
                    {
                        const auto& intf = interfaces_[inti].interface();
                        const auto* ppp =
                            isA<const processorLduInterface>(intf);
                        if (ppp)
                        {
                            if
                            (
                                colours[Pstream::myProcNo()] == colouri
                             && colours[ppp->neighbProcNo()] > colouri
                            )
                            {
                                // Neighbour not yet updated
                                meFromHigherNbrs.append(inti);
                            }
                            else if
                            (
                                colours[Pstream::myProcNo()] == colouri
                             && colours[ppp->neighbProcNo()] < colouri
                            )
                            {
                                // Neighbour not yet updated
                                meFromLowerNbrs.append(inti);
                            }
                            if
                            (
                                colours[Pstream::myProcNo()] > colouri
                             && colours[ppp->neighbProcNo()] == colouri
                            )
                            {
                                // Neighbour just been solved for. I am not
                                // yet solved for
                                higherFromMeNbrs.append(inti);
                            }
                        }
                    }
                }
                //DebugVar(meFromHigherNbrs);
                //DebugVar(higherFromMeNbrs);


                // Include not-solved for neighbours
                scalarField rD_;
                calcReciprocalD(rD_, matrix(), meFromHigherNbrs, coeffs, cmpt);

                //DebugVar(rD_);



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

                // Initialise 'internal' cells
                if (colours[Pstream::myProcNo()] == colouri)
                {
                    for (label cell=0; cell<nCells; cell++)
                    {
                        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
                    }
                }

                // Do 'halo' contributions from higher neighbours
                {
                    forAll(interfaces_, inti)
                    {
                        if (interfaces_.set(inti))
                        {
                            coeffs[inti] = Zero;
                        }
                    }
                    for (const label inti : meFromHigherNbrs)
                    {
                        const auto& intf = interfaces_[inti].interface();
                        const auto& faceCells = intf.faceCells();
                        const auto& bc = interfaceBouCoeffs_[inti];
                        scalarField& coeff = coeffs[inti];

                        forAll(bc, i)
                        {
                            // bouCoeffs are negative compared to
                            // upperPtr
                            coeff[i] = -bc[i]*rD_[faceCells[i]];
                        }
                    }

                    const label startRequest = Pstream::nRequests();
                    matrix().initMatrixInterfaces
                    (
                        true,  // subtract contribution
                        coeffs,
                        interfaces_,
                        wA,
                        wA,
                        cmpt                
                    );
                    matrix().updateMatrixInterfaces
                    (
                        true,  // subtract contribution
                        coeffs,
                        interfaces_,
                        wA,
                        wA,
                        cmpt,
                        startRequest              
                    );
                }

                if (colours[Pstream::myProcNo()] == colouri)
                {
                    // Do 'internal' faces
                    for (label face=0; face<nFaces; face++)
                    {
                        wAPtr[uPtr[face]] -=
                            rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
                    }
                }

                // Do 'halo' contributions from lower neighbours
                if (false)
                {
                    forAll(interfaces_, inti)
                    {
                        if (interfaces_.set(inti))
                        {
                            coeffs[inti] = Zero;
                        }
                    }
                    for (const label inti : meFromLowerNbrs)
                    {
                        const auto& intf = interfaces_[inti].interface();
                        const auto& faceCells = intf.faceCells();
                        const auto& bc = interfaceBouCoeffs_[inti];
                        scalarField& coeff = coeffs[inti];

                        forAll(bc, i)
                        {
                            // bouCoeffs are negative compared to
                            // upperPtr
                            coeff[i] = -bc[i]*rD_[faceCells[i]];
                        }
                    }

                    const label startRequest = Pstream::nRequests();
                    matrix().initMatrixInterfaces
                    (
                        true,  // subtract contribution
                        coeffs,
                        interfaces_,
                        wA,
                        wA,
                        cmpt                
                    );
                    matrix().updateMatrixInterfaces
                    (
                        true,  // subtract contribution
                        coeffs,
                        interfaces_,
                        wA,
                        wA,
                        cmpt,
                        startRequest              
                    );
                }

                if (colours[Pstream::myProcNo()] == colouri)
                {
                    for (label face=nFacesM1; face>=0; face--)
                    {
                        wAPtr[lPtr[face]] -=
                            rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
                    }
                }

//                // Send back to origating processors
//                forAll(interfaces_, inti)
//                {
//                    if (interfaces_.set(inti))
//                    {
//                        coeffs[inti] = Zero;
//                    }
//                }
//                for (const label inti : higherFromMeNbrs)
//                {
//                    const auto& intf = interfaces_[inti].interface();
//                    const auto& faceCells = intf.faceCells();
//                    const auto& bc = interfaceBouCoeffs_[inti];
//                    scalarField& coeff = coeffs[inti];
//
//                    forAll(bc, i)
//                    {
//                        // bouCoeffs are negative compared
//                        // to upperPtr
//                        coeff[i] = -bc[i]*rD_[faceCells[i]];
//                    }
//                }
//
//                const label oldRequest = Pstream::nRequests();
//                matrix().initMatrixInterfaces
//                (
//                    true,  // subtract contribution
//                    coeffs,
//                    interfaces_,
//                    wA,
//                    wA,
//                    cmpt                
//                );
//                matrix().updateMatrixInterfaces
//                (
//                    true,  // subtract contribution
//                    coeffs,
//                    interfaces_,
//                    wA,
//                    wA,
//                    cmpt,
//                    oldRequest              
//                );
            }

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
