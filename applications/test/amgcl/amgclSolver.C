/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "amgcl/amg.hpp"
#include "amgcl/make_solver.hpp"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/adapter/crs_tuple.hpp"
//#include "amgcl/coarsening/ruge_stuben.hpp"
//#include "amgcl/coarsening/pointwise_aggregates.hpp"
#include "amgcl/coarsening/aggregation.hpp"
#include "amgcl/coarsening/smoothed_aggregation.hpp"
//#include "amgcl/relaxation/damped_jacobi.hpp"
#include "amgcl/relaxation/spai0.hpp"
#include "amgcl/solver/bicgstab.hpp"
//#include "amgcl/solver/gmres.hpp"

#include "amgclSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(amgclSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<amgclSolver>
        addamgclSolverSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<amgclSolver>
        addamgclSolverAsymMatrixConstructorToTable_;
}


typedef amgcl::backend::builtin<double> Backend;

typedef amgcl::make_solver
<
    // Use AMG as preconditioner:
    amgcl::amg
    <
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
        >,
    // And BiCGStab as iterative solver:
    amgcl::solver::bicgstab<Backend>
> Solver;


void Foam::amgclSolver::crs
(
    lduMatrix& matrix,
    std::vector<double>& val,
    std::vector<int>& col,
    std::vector<int>& ptr
    //std::vector<double>& rhs
)
{
    const lduAddressing& lduAddr = matrix.mesh().lduAddr();
    const labelUList& uPtr = lduAddr.upperAddr();
    const labelUList& lPtr = lduAddr.lowerAddr();
    const label n = lduAddr.size();
    const scalarField& diag = matrix.diag();
    const scalarField& upperPtr = matrix.upper();

    const scalarField& lowerPtr =
    (
        matrix.hasLower()
      ? matrix.lower()
      : matrix.upper()
    );

    // Count number of neighbours
    labelList nNbrs(n, 0);
    forAll(uPtr, facei)
    {
        nNbrs[uPtr[facei]]++;
        nNbrs[lPtr[facei]]++;
    }

    // Points to the start of each row in the above arrays
    ptr.resize(n+1);

    ptr[0] = 0;
    for (label celli = 0; celli < n; celli++)
    {
        // Space for diagonal and cell neighbours
        ptr[celli+1] = ptr[celli]+1+nNbrs[celli];
    }
    // Values of nonzero entries
    val.resize(ptr[ptr.size()-1]);
    // Column numbers of nonzero entries
    col.resize(val.size());

    std::vector<int> offsets(ptr);
    // Store diagonal first (not needed per se but easier when adding
    // boundary contributions
    for (label celli = 0; celli < n; celli++)
    {
        int& index = offsets[celli];
        col[index] = celli;
        val[index++] = diag[celli];
    }
    // Store connections to lower numbered cells
    forAll(uPtr, facei)
    {
        int& index = offsets[uPtr[facei]];
        col[index] = lPtr[facei];
        val[index++] = lowerPtr[facei];
    }

    // Store connections to higher numbered cells
    forAll(lPtr, facei)
    {
        int& index = offsets[lPtr[facei]];
        col[index] = uPtr[facei];
        val[index++] = upperPtr[facei];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::amgclSolver::amgclSolver
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
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::amgclSolver::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    scalar normFactor = 0;

    {
        scalarField Apsi(psi.size());
        scalarField temp(psi.size());

        // Calculate A.psi
        matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

        // Calculate normalisation factor
        normFactor = this->normFactor(psi, source, Apsi, temp);

        // Calculate residual magnitude
        solverPerf.initialResidual() = gSumMag
        (
            (source - Apsi)(),
            matrix().mesh().comm()
        )/normFactor;
        solverPerf.finalResidual() = solverPerf.initialResidual();
    }

    if (lduMatrix::debug >= 2)
    {
        Info.masterStream(matrix().mesh().comm())
            << "   Normalisation factor = " << normFactor << endl;
    }


    // Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
        // Convert to CRS format
        std::vector<double> val;
        std::vector<int> col;
        std::vector<int> ptr;
        crs(const_cast<lduMatrix&>(matrix()), val, col, ptr);

        // Do rhs
        const lduAddressing& lduAddr = matrix().mesh().lduAddr();
        const label n = lduAddr.size();

        std::vector<double> rhs(n);
        forAll(source, celli)
        {
            rhs[celli] = source[celli];
        }
        forAll(interfaceBouCoeffs(), patchi)
        {
            const scalarField& iCoeffs = interfaceIntCoeffs()[patchi];
            const scalarField& bCoeffs = interfaceBouCoeffs()[patchi];
            const labelUList& fc = lduAddr.patchAddr(patchi);

            forAll(fc, i)
            {
                label celli = fc[i];
                label index = ptr[celli];
                if (col[index] != celli)
                {
                    FatalErrorInFunction
                        << "Problem: diagonal not first" << abort(FatalError);
                }

                val[index] += iCoeffs[i];
                rhs[celli] += bCoeffs[i];
            }
        }


        // Construct the iterative solver. It needs size of the system to
        // preallocate the required temporary structures:
        Solver::params prm;
        prm.solver.tol = tolerance_;    ///normFactor;
        Solver solve(boost::tie(n, ptr, col, val), prm);

        // Solve the system. Returns number of iterations made and the
        // achieved residual.
        std::vector<double> solution(n);
        forAll(psi, celli)
        {
            solution[celli] = psi[celli];
        }
        label iters;
        scalar resid;
        boost::tie(iters, resid) = solve(rhs, solution);
        //Info << iters << "   " << resid << endl;

        forAll(psi, celli)
        {
            psi[celli] = solution[celli];
        }


        // Calculate the residual to check convergence
        solverPerf.nIterations() = iters;
        solverPerf.finalResidual() = gSumMag
        (
            matrix_.residual
            (
                psi,
                source,
                interfaceBouCoeffs_,
                interfaces_,
                cmpt
            )(),
            matrix().mesh().comm()
        )/normFactor;
    }

    return solverPerf;
}


// ************************************************************************* //
