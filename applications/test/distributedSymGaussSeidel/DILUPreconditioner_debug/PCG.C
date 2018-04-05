/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "PCG.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PCG_debug, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PCG_debug>
        addPCG_debugSymMatrixConstructorToTable_;
}


void Foam::PCG_debug::updateMatrixInterfaces
(
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const labelUList& selectedInterfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt
) const
{
    forAll(selectedInterfaces, i)
    {
        label interfacei = selectedInterfaces[i];

        interfaces[interfacei].initInterfaceMatrixUpdate
        (
            result,
            psiif,
            coupleCoeffs[interfacei],
            cmpt,
            Pstream::defaultCommsType
        );
    }

    // Block for all requests and remove storage
    UPstream::waitRequests();

    forAll(selectedInterfaces, i)
    {
        label interfacei = selectedInterfaces[i];
        interfaces[interfacei].updateInterfaceMatrix
        (
            result,
            psiif,
            coupleCoeffs[interfacei],
            cmpt,
            Pstream::defaultCommsType
        );

        Pout<< "   interface:" << interfacei
            << " adding contributions with weigts:" << coupleCoeffs[interfacei]
            << " to cells:" << interfaces[interfacei].interface().faceCells()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PCG_debug::PCG_debug
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

Foam::solverPerformance Foam::PCG_debug::solve
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


    // Pre-calculate 'left' and 'right' processor interfaces
    DynamicList<label> lowerInterfaces(interfaces_.size());
    DynamicList<label> upperInterfaces(interfaces_.size());
    forAll(interfaces_, inti)
    {
        if
        (
            interfaces_.set(inti)
         && isA<processorLduInterface>(interfaces_[inti].interface())
        )
        {
            const processorLduInterface& procInt =
            refCast<const processorLduInterface>
            (
                interfaces_[inti].interface()
            );

            Pout<< "    interface:" << inti
                << " type:" << interfaces_[inti].interface().type()
                << " myProcNo:" << procInt.myProcNo()
                << " neighbProcNo:" << procInt.neighbProcNo()
                << endl;

            if (procInt.neighbProcNo() > procInt.myProcNo())
            {
                upperInterfaces.append(inti);
            }
            else
            {
                lowerInterfaces.append(inti);
            }
        }
    }
    Pout<< "lowerInterfaces:" << lowerInterfaces << endl;
    Pout<< "upperInterfaces:" << upperInterfaces << endl;

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


    const label nFaces = matrix().upper().size();
    const label* const __restrict__ uPtr =
        matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix().lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix().upper().begin();

DebugVar(matrix().upper());
DebugVar(interfaceBouCoeffs_);


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
        scalarField rD_(matrix().diag());

        {
            // Subtract coupled contributions
            scalarField invRd(matrix().diag());
            forAll(invRd, cell)
            {
                invRd[cell] = 1.0/invRd[cell];
            }
            FieldField<Field, scalar> coeffs
            (
                1.0*interfaceBouCoeffs_*interfaceBouCoeffs_
            );
            updateMatrixInterfaces
            (
                coeffs,
                interfaces_,
                lowerInterfaces,
                invRd,
                rD_,
                cmpt
            );

            forAll(rD_, cell)
            {
                Pout<< "diag:" << matrix().diag()[cell]
                    << " rD_:" << rD_[cell] << endl;
            }
        }
        {
            for (label face=0; face<nFaces; face++)
            {
                Pout<< "For face:" << face
                    << " adapting diagonal for cell:" << uPtr[face]
                    << " contributions from cell:" << lPtr[face]
                    << " from " << rD_[uPtr[face]];

                rD_[uPtr[face]] -=
                    upperPtr[face]*upperPtr[face]/rD_[lPtr[face]];

                Pout<< " to " << rD_[uPtr[face]] << endl;
            }


            // Calculate the reciprocal of the preconditioned diagonal
            const label nCells = rD_.size();

            for (label cell=0; cell<nCells; cell++)
            {
                rD_[cell] = 1.0/rD_[cell];
            }
        }


        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
            //preconPtr->precondition(wA, rA, cmpt);
            {
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

                for (label cell=0; cell<nCells; cell++)
                {
                    wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
                }


                // Calculate coupled contributions
                //if (false)
                {
                    FieldField<Field, scalar> coeffs
                    (
                        1.0*interfaceBouCoeffs_
                    );
                    forAll(lowerInterfaces, i)
                    {
                        label inti = lowerInterfaces[i];
                        coeffs[inti] *=
                            scalarField
                            (
                                rD_,
                                interfaces_[inti].interface().faceCells()
                            );
                    }
                    updateMatrixInterfaces
                    (
                        coeffs,
                        interfaces_,
                        lowerInterfaces,
                        wA,
                        wA,
                        cmpt
                    );
                }


                for (label face=0; face<nFaces; face++)
                {
                    Pout<< "For face:" << face
                        << " adapting cell:" << uPtr[face]
                        << " for contributions from cell:" << lPtr[face]
                        << " from:" << wAPtr[uPtr[face]];

                    wAPtr[uPtr[face]] -=
                        rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];

                    Pout<< " to " << wAPtr[uPtr[face]] << endl;
                }

                {
                    FieldField<Field, scalar> coeffs
                    (
                        1.0*interfaceBouCoeffs_
                    );
                    forAll(upperInterfaces, i)
                    {
                        label inti = upperInterfaces[i];
                        coeffs[inti] *=
                            scalarField
                            (
                                rD_,
                                interfaces_[inti].interface().faceCells()
                            );
                    }
                    updateMatrixInterfaces
                    (
                        coeffs,
                        interfaces_,
                        upperInterfaces,
                        wA,
                        wA,
                        cmpt
                    );
                }


                for (label face=nFacesM1; face>=0; face--)
                {
                    Pout<< "For face:" << face
                        << " adapting cell:" << lPtr[face]
                        << " for contributions from cell:" << uPtr[face]
                        << " from:" << wAPtr[lPtr[face]];

                    wAPtr[lPtr[face]] -=
                        rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];

                    Pout<< " to " << wAPtr[lPtr[face]] << endl;
                }
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
