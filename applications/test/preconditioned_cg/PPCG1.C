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

#include "PPCG1.H"
//#include <mpi.h>
//
//#if defined(WM_SP)
//    #define MPI_SCALAR MPI_FLOAT
//#elif defined(WM_DP)
//    #define MPI_SCALAR MPI_DOUBLE
//#elif defined(WM_LP)
//    #define MPI_SCALAR MPI_LONG_DOUBLE
//#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PPCG1, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PPCG1>
        addPPCG1SymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PPCG1::calcDirections
(
    FixedList<scalar, 3>& globalSum,
    const scalarField& r,
    const scalarField& u,
    const scalarField& w,
    label& outstandingRequest,
    const label comm
) const
{
    const label nCells = r.size();

    globalSum = 0.0;
    for (label cell=0; cell<nCells; cell++)
    {
        globalSum[0] += w[cell]*u[cell];
        globalSum[1] += r[cell]*u[cell];

        //- Convergence based on preconditioned residual
        //globalSum[2] += mag(u[cell]);
        //- Convergence based on non-preconditioned residual
        globalSum[2] += mag(r[cell]);
    }

    if (Pstream::parRun())
    {
//        const int err = MPI_Iallreduce
//        (
//            MPI_IN_PLACE,       //globalSum.cbegin(),
//            globalSum.begin(),
//            globalSum.size(),   //MPICount,
//            MPI_SCALAR,         //MPIType,
//            MPI_SUM,            //MPIOp,
//            MPI_COMM_WORLD,     //TBD. comm,
//            &outstandingRequest
//        );
//        if (err)
//        {
//            FatalErrorInFunction<< "Failed MPI_Iallreduce for "
//                << globalSum << exit(FatalError);
//        }
        Foam::reduce
        (
            globalSum.begin(),
            globalSum.size(),
            Pstream::msgType(),
            comm,
            outstandingRequest
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PPCG1::PPCG1
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

Foam::solverPerformance Foam::PPCG1::solve
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
    scalarField wA(nCells);

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    scalarField r(source - wA);

    // --- Calculate normalisation factor
    scalarField p(nCells);
    const scalar normFactor = this->normFactor(psi, source, wA, p);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

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


    // State
    scalarField s(nCells);
    scalarField q(nCells);
    scalarField z(nCells);


    // --- Start global reductions for inner products
    FixedList<scalar, 3> globalSum;
    label outstandingRequest;
    calcDirections(globalSum, r, u, w, outstandingRequest, comm);

    scalar alpha = 0.0;


    // --- Precondition residual
    scalarField m(nCells);
    preconPtr->precondition(m, w, cmpt);

    // --- Calculate A*m
    scalarField n(nCells);
    matrix_.Amul(n, m, interfaceBouCoeffs_, interfaces_, cmpt);

    scalar gamma = 0.0;

    // --- Solver iteration
    for
    (
        solverPerf.nIterations() = 0;
        solverPerf.nIterations() < maxIter_;
        solverPerf.nIterations()++
    )
    {
        // Make sure gamma,delta are available
        if (Pstream::parRun())
        {
            Pstream::waitRequest(outstandingRequest);
        }

        const scalar delta = globalSum[0];
        const scalar gammaOld = gamma;
        gamma = globalSum[1];

        solverPerf.finalResidual() = globalSum[2]/normFactor;
        if (solverPerf.nIterations() == 0)
        {
            solverPerf.initialResidual() = solverPerf.finalResidual();
        }

        // Check convergence (bypass if not enough iterations yet)
        if
        (
            (minIter_ <= 0 || solverPerf.nIterations() >= minIter_)
         && solverPerf.checkConvergence(tolerance_, relTol_)
        )
        {
            break;
        }


        if (solverPerf.nIterations() == 0)
        {
            alpha = gamma/delta;

            z = n;
            q = m;
            s = w;
            p = u;
        }
        else
        {
            const scalar beta = gamma/gammaOld;
            alpha = gamma/(delta-beta*gamma/alpha);

            for (label cell=0; cell<nCells; cell++)
            {
                z[cell] = n[cell] + beta*z[cell];
                q[cell] = m[cell] + beta*q[cell];
                s[cell] = w[cell] + beta*s[cell];
                p[cell] = u[cell] + beta*p[cell];
            }
        }

        for (label cell=0; cell<nCells; cell++)
        {
            psi[cell] += alpha*p[cell];
            r[cell] -= alpha*s[cell];
            u[cell] -= alpha*q[cell];
            w[cell] -= alpha*z[cell];
        }

        // --- Start global reductions for inner products
        calcDirections(globalSum, r, u, w, outstandingRequest, comm);

        // --- Precondition residual
        preconPtr->precondition(m, w, cmpt);

        // --- Calculate A*m
        matrix_.Amul(n, m, interfaceBouCoeffs_, interfaces_, cmpt);
    }

    return solverPerf;
}


// ************************************************************************* //
