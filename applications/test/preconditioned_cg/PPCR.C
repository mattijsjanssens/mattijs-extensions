/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "PPCR.H"
#include <mpi.h>

#if defined(WM_SP)
    #define MPI_SCALAR MPI_FLOAT
#elif defined(WM_DP)
    #define MPI_SCALAR MPI_DOUBLE
#elif defined(WM_LP)
    #define MPI_SCALAR MPI_LONG_DOUBLE
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PPCR, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PPCR>
        addPPCRSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PPCR::PPCR
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

Foam::solverPerformance Foam::PPCR::solve
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

    const label comm = matrix().mesh().comm();
    const label nCells = psi.size();
    scalarField p(nCells);
    scalarField wA(nCells);

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    scalarField r(source - wA);

    // --- Calculate normalisation factor
    scalar normFactor = this->normFactor(psi, source, wA, p);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() = gSumMag(r, comm)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
        // --- Select and construct the preconditioner
        autoPtr<lduMatrix::preconditioner> preconPtr =
        lduMatrix::preconditioner::New
        (
            *this,
            controlDict_
        );


        // --- Precondition residual (= u0)
        scalarField u(nCells);
        preconPtr->precondition(u, r, cmpt);

        // --- Calculate A*u
        scalarField w(nCells);
        matrix_.Amul(w, u, interfaceBouCoeffs_, interfaces_, cmpt);


        scalarField m(nCells);
        scalarField n(nCells);
        scalarField q(nCells);
        scalarField z(nCells);

        scalar gamma = 0.0;
        scalar gammaOld;
        scalar alpha = 0.0;
        scalar alphaOld;

        scalarField localGammaDelta(2);
        scalarField gammaDelta(2);
        MPI_Request outstandingRequest;

        // --- Solver iteration
        do
        {
            // --- Precondition residual
            preconPtr->precondition(m, w, cmpt);

            // --- Update search directions:
            gammaOld = gamma;

            //gamma = gSumProd(w, u, comm);
            localGammaDelta[0] = sumProd(w, u);
            localGammaDelta[1] = sumProd(m, w);

            if (Pstream::parRun())
            {
                const int err = MPI_Iallreduce
                (
                    localGammaDelta.cbegin(),
                    gammaDelta.begin(),
                    int(2),     //MPICount,
                    MPI_SCALAR, //MPIType,
                    MPI_SUM,    //MPIOp,
                    MPI_COMM_WORLD,          //TBD. comm,
                    &outstandingRequest
                );
                if (err)
                {
                    FatalErrorInFunction<< "Failed MPI_Iallreduce for "
                        << localGammaDelta << exit(FatalError);
                }
            }
            else
            {
                gammaDelta = localGammaDelta;
            }

            matrix_.Amul(n, m, interfaceBouCoeffs_, interfaces_, cmpt);

            // Make sure gamma,delta are available
            if (Pstream::parRun())
            {
                if (MPI_Wait(&outstandingRequest, MPI_STATUS_IGNORE))
                {
                    FatalErrorInFunction<< "Failed waiting for"
                        << " MPI_Iallreduce request" << exit(FatalError);
                }
            }
            gamma = gammaDelta[0];
            const scalar delta = gammaDelta[1];

            alphaOld = alpha;
            if (solverPerf.nIterations() == 0)
            {
                alpha = gamma/delta;
                z = n;
                q = m;
                p = u;
            }
            else
            {
                const scalar beta = gamma/gammaOld;
                alpha = gamma/(delta-beta*gamma/alphaOld);

                z = n + beta*z;
                q = m + beta*q;
                p = u + beta*p;
            }

            for (label cell=0; cell<nCells; cell++)
            {
                psi[cell] += alpha*p[cell];
                u[cell] -= alpha*q[cell];
                w[cell] -= alpha*z[cell];
            }

            matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);
            r = source - wA;
            solverPerf.finalResidual() =
                gSumMag(r, comm)
               /normFactor;

        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    return solverPerf;
}


// ************************************************************************* //
