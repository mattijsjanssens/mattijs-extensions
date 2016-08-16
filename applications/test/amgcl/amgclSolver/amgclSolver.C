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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/scope_exit.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/algorithm.hpp>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/multi_array.hpp>

#include <amgcl/make_solver.hpp>
#include <amgcl/runtime.hpp>
#include <amgcl/mpi/direct_solver.hpp>
#include <amgcl/mpi/subdomain_deflation.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/profiler.hpp>

#include "amgclSolver.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(amgclSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<amgclSolver>
        addamgclSolverSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<amgclSolver>
        addamgclSolverAsymMatrixConstructorToTable_;
}

struct constant_deflation {
    size_t dim() const { return 1; }
    double operator()(ptrdiff_t i, unsigned j) const { return 1.0; }
};


void Foam::amgclSolver::crs
(
    const globalIndex& globalNumbering,
    lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const scalarField& source,
    std::vector<double>& val,
    std::vector<int>& col,
    std::vector<int>& ptr,
    std::vector<double>& rhs
)
{
    const lduAddressing& lduAddr = matrix.mesh().lduAddr();
    const lduInterfacePtrsList interfaces(matrix.mesh().interfaces());
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
    forAll(interfaces, patchi)
    {
        // Is coupled interface with neighbour?
        if (interfaces.set(patchi))
        {
            const labelUList& fc = lduAddr.patchAddr(patchi);
            forAll(fc, i)
            {
                nNbrs[fc[i]]++;
            }
        }
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
    // boundary contributions)
    for (label celli = 0; celli < n; celli++)
    {
        int& index = offsets[celli];
        col[index] = globalNumbering.toGlobal(celli);
        val[index++] = diag[celli];
    }
    // Store connections to lower numbered cells
    forAll(uPtr, facei)
    {
        int& index = offsets[uPtr[facei]];
        col[index] = globalNumbering.toGlobal(lPtr[facei]);
        val[index++] = lowerPtr[facei];
    }

    // Store connections to higher numbered cells
    forAll(lPtr, facei)
    {
        int& index = offsets[lPtr[facei]];
        col[index] = globalNumbering.toGlobal(uPtr[facei]);
        val[index++] = upperPtr[facei];
    }

    // Store connections to neighbouring processors
    {
        labelList globalCells(n);
        forAll(globalCells, celli)
        {
            globalCells[celli] = globalNumbering.toGlobal(celli);
        }

        // Initialise transfer of global cells
        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi))
            {
                interfaces[patchi].initInternalFieldTransfer
                (
                    Pstream::nonBlocking,
                    globalCells
                );
            }
        }

        if (Pstream::parRun())
        {
            Pstream::waitRequests();
        }

        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi))
            {
                labelField nbrCells
                (
                    interfaces[patchi].internalFieldTransfer
                    (
                        Pstream::nonBlocking,
                        globalCells
                    )
                );

                const labelUList& fc = lduAddr.patchAddr(patchi);
                const scalarField& bCoeffs = interfaceBouCoeffs[patchi];
                forAll(fc, i)
                {
                    label celli = fc[i];
                    int& index = offsets[celli];
                    col[index] = nbrCells[i];

                    // Note: opposite sign since we're using this sides'
                    //       coefficients (from discretising our face; not the
                    //       neighbouring, reversed face)
                    val[index++] = -bCoeffs[i];
                }
            }
        }
    }

    // Do rhs
    rhs.resize(n);
    forAll(source, celli)
    {
        rhs[celli] = source[celli];
    }
    forAll(interfaceBouCoeffs, patchi)
    {
        forAll(interfaces, patchi)
        {
            if (!interfaces.set(patchi))
            {
                const scalarField& iCoeffs = interfaceIntCoeffs[patchi];
                const scalarField& bCoeffs = interfaceBouCoeffs[patchi];
                const labelUList& fc = lduAddr.patchAddr(patchi);

                forAll(fc, i)
                {
                    label celli = fc[i];
                    label index = ptr[celli];
                    if (col[index] != globalNumbering.toGlobal(celli))
                    {
                        FatalErrorInFunction
                            << "Problem: diagonal not first"
                            << abort(FatalError);
                    }

                    val[index] += iCoeffs[i];
                    rhs[celli] += bCoeffs[i];
                }
            }
        }
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
        globalIndex globalNumbering(matrix().mesh().lduAddr().size());

        // Convert to CRS format
        std::vector<double> val;
        std::vector<int> col;
        std::vector<int> ptr;
        std::vector<double> rhs;
        crs
        (
            globalNumbering,
            const_cast<lduMatrix&>(matrix()),
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            source,
            val,
            col,
            ptr,
            rhs
        );


        std::vector<double> solution(rhs.size());
        forAll(psi, celli)
        {
            solution[celli] = psi[celli];
        }


        typedef amgcl::backend::builtin<double> Backend;


        ptrdiff_t n(globalNumbering.localSize());

        size_t iters;
        double resid;
        if (Pstream::parRun())
        {
            //amgcl::runtime::coarsening::type coarsening =
            //    amgcl::runtime::coarsening::smoothed_aggregation;
            //amgcl::runtime::relaxation::type relaxation =
            //    amgcl::runtime::relaxation::spai0;
            //amgcl::runtime::solver::type iterative_solver =
            //    amgcl::runtime::solver::bicgstabl;
            //amgcl::runtime::mpi::dsolver::type direct_solver =
            //    amgcl::runtime::mpi::dsolver::skyline_lu;
            //
            //
            //boost::property_tree::ptree prm;
            //prm.put("isolver.type", iterative_solver);
            //prm.put("dsolver.type", direct_solver);
            //prm.put("solver.tol", tolerance_);
            //prm.put("isolver.tol", tolerance_);
            //
            //boost::function<double(ptrdiff_t,unsigned)> dv;
            //unsigned ndv;
            //
            //// Use constant deflation for now
            //{
            //    dv = constant_deflation();
            //    ndv = 1;
            //}
            //prm.put("num_def_vec", ndv);
            //prm.put("def_vec", &dv);
            //
            //prm.put("local.coarsening.type", coarsening);
            //prm.put("local.relax.type", relaxation);
            //    typedef
            //        amgcl::mpi::subdomain_deflation
            //        <
            //            amgcl::runtime::amg<Backend>,
            //            amgcl::runtime::iterative_solver,
            //            amgcl::runtime::mpi::direct_solver<double>
            //        > Solver;

            typedef amgcl::mpi::subdomain_deflation
            <
                // Use AMG as preconditioner:
                amgcl::amg
                <
                    Backend,
                    amgcl::coarsening::smoothed_aggregation,
                    amgcl::relaxation::spai0
                >,
                // Iterative solver
                amgcl::solver::bicgstab,
                // Direct solver
                amgcl::mpi::skyline_lu<double>
            > Solver;

            Solver::params prm;
            prm.isolver.tol = tolerance_;
            prm.num_def_vec = 1;
            prm.def_vec = constant_deflation();

            amgcl::mpi::communicator world(MPI_COMM_WORLD);

            Solver solve(world, boost::tie(n, ptr, col, val), prm);

            boost::tie(iters, resid) = solve(rhs, solution);
        }
        else
        {
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

            Solver::params prm;
            prm.solver.tol = tolerance_;

            Solver solve(boost::tie(n, ptr, col, val), prm);

            boost::tie(iters, resid) = solve(rhs, solution);
        }


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
