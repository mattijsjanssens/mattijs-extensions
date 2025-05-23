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

#include "PPCG_Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PPCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PPCG>
        addPPCGSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PPCG::gSumMagProd
(
    FixedList<scalar, 3>& globalSum,
    const scalarField& a,
    const scalarField& b,
    const scalarField& c,
    const scalarField& sumMag,
    MPI_Request& outstandingRequest
) const
{
    const label nCells = a.size();

    globalSum = 0.0;
    for (label cell=0; cell<nCells; cell++)
    {
        globalSum[0] += a[cell]*b[cell];    // sumProd(a, b)
        globalSum[1] += a[cell]*c[cell];    // sumProd(a, c)
        globalSum[2] += mag(sumMag[cell]);
    }

    if (Pstream::parRun())
    {
        const int err = MPI_Iallreduce
        (
            MPI_IN_PLACE,       //globalSum.cbegin(),
            globalSum.begin(),
            globalSum.size(),   //MPICount,
            MPI_SCALAR,         //MPIType,
            MPI_SUM,            //MPIOp,
            MPI_COMM_WORLD,     //TBD. comm,
            &outstandingRequest
        );
        if (err)
        {
            FatalErrorInFunction<< "Failed MPI_Iallreduce for "
                << globalSum << exit(FatalError);
        }
    }
}


Foam::solverPerformance Foam::PPCG::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const bool cgMode
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    const label nCells = psi.size();
    scalarField w(nCells);

    // --- Calculate A.psi
    matrix_.Amul(w, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    scalarField r(source - w);

    // --- Calculate normalisation factor
    scalarField p(nCells);
    const scalar normFactor = this->normFactor(psi, source, w, p);

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

    // --- Calculate A*u - reuse w
    matrix_.Amul(w, u, interfaceBouCoeffs_, interfaces_, cmpt);


    // State
    scalarField s(nCells);
    scalarField q(nCells);
    scalarField z(nCells);

    scalarField m(nCells);

    FixedList<scalar, 3> globalSum;
    MPI_Request outstandingRequest;
    if (cgMode)
    {
        // --- Start global reductions for inner products
        gSumMagProd(globalSum, u, r, w, r, outstandingRequest);

        // --- Precondition residual
        preconPtr->precondition(m, w, cmpt);
    }
    else
    {
        // --- Precondition residual
        preconPtr->precondition(m, w, cmpt);

        // --- Start global reductions for inner products
        gSumMagProd(globalSum, w, u, m, r, outstandingRequest);
    }

    // --- Calculate A*m
    scalarField n(nCells);
    matrix_.Amul(n, m, interfaceBouCoeffs_, interfaces_, cmpt);

    scalar alpha = 0.0;
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
            if (MPI_Wait(&outstandingRequest, MPI_STATUS_IGNORE))
            {
                FatalErrorInFunction<< "Failed waiting for"
                    << " MPI_Iallreduce request" << exit(FatalError);
            }
        }

        const scalar gammaOld = gamma;
        gamma = globalSum[0];
        const scalar delta = globalSum[1];

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

        if (cgMode)
        {
            // --- Start global reductions for inner products
            gSumMagProd(globalSum, u, r, w, r, outstandingRequest);

            // --- Precondition residual
            preconPtr->precondition(m, w, cmpt);
        }
        else
        {
            // --- Precondition residual
            preconPtr->precondition(m, w, cmpt);

            // --- Start global reductions for inner products
            gSumMagProd(globalSum, w, u, m, r, outstandingRequest);
        }

        // --- Calculate A*m
        matrix_.Amul(n, m, interfaceBouCoeffs_, interfaces_, cmpt);
    }

    return solverPerf;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PPCG::PPCG
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

Foam::solverPerformance Foam::PPCG::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    return solve(psi, source, cmpt, true);
}


// ************************************************************************* //
