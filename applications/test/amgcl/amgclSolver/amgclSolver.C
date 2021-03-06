/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include <amgcl/amg.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/runtime.hpp>
//#include <amgcl/mpi/direct_solver.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/mpi/subdomain_deflation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/skyline_lu.hpp>
#include <amgcl/profiler.hpp>

#include "amgclSolver.H"
#include "globalIndex.H"
#include "cpuTime.H"
#include "OFstream.H"

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
    const scalarField& source,
    std::vector<double>& val,
    std::vector<int>& col,
    std::vector<double>& rhs
) const
{
    const lduAddressing& lduAddr = matrix_.mesh().lduAddr();
    const lduInterfacePtrsList interfaces(matrix_.mesh().interfaces());
    const labelUList& uPtr = lduAddr.upperAddr();
    const labelUList& lPtr = lduAddr.lowerAddr();
    const label n = lduAddr.size();
    const scalarField& diag = matrix_.diag();
    const scalarField& upperPtr = matrix_.upper();

    const scalarField& lowerPtr =
    (
        matrix_.hasLower()
      ? matrix_.lower()
      : matrix_.upper()
    );

    // Values of nonzero entries
    val.resize(ptr_[ptr_.size()-1]);
    // Column numbers of nonzero entries
    col.resize(val.size());

    std::vector<int> offsets(ptr_);
    // Store diagonal first (not needed per se but easier when adding
    // boundary contributions)
    for (label celli = 0; celli < n; celli++)
    {
        int& index = offsets[celli];
        col[index] = globalNumbering_.toGlobal(celli);
        val[index++] = diag[celli];
    }
    // Store connections to lower numbered cells
    forAll(uPtr, facei)
    {
        int& index = offsets[uPtr[facei]];
        col[index] = globalNumbering_.toGlobal(lPtr[facei]);
        val[index++] = lowerPtr[facei];
    }

    // Store connections to higher numbered cells
    forAll(lPtr, facei)
    {
        int& index = offsets[lPtr[facei]];
        col[index] = globalNumbering_.toGlobal(uPtr[facei]);
        val[index++] = upperPtr[facei];
    }

    // Store connections to neighbouring processors
    {
        labelList globalCells(n);
        forAll(globalCells, celli)
        {
            globalCells[celli] = globalNumbering_.toGlobal(celli);
        }

        // Initialise transfer of global cells
        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi))
            {
                interfaces[patchi].initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
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
                        Pstream::commsTypes::nonBlocking,
                        globalCells
                    )
                );

                const labelUList& fc = lduAddr.patchAddr(patchi);
                const scalarField& bCoeffs = interfaceBouCoeffs_[patchi];
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
    ),
    globalNumbering_(matrix.mesh().lduAddr().size())
{
    readControls();

    const lduAddressing& lduAddr = matrix.mesh().lduAddr();
    const labelUList& uPtr = lduAddr.upperAddr();
    const labelUList& lPtr = lduAddr.lowerAddr();
    const label n = lduAddr.size();

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
    ptr_.resize(n+1);

    ptr_[0] = 0;
    for (label celli = 0; celli < n; celli++)
    {
        // Space for diagonal and cell neighbours
        ptr_[celli+1] = ptr_[celli]+1+nNbrs[celli];
    }
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

        cpuTime solverTime;

        scalar tol = max
        (
            relTol_*solverPerf.initialResidual(),
            tolerance_
        );


        // Convert to CRS format
        std::vector<double> val;
        std::vector<int> col;
        std::vector<double> rhs;
        crs(source, val, col, rhs);

        if (debug)
        {
            Pout<< "Constructed addressing in " << solverTime.cpuTimeIncrement()
                << endl;
        }

        if (debug & 2)
        {
            // Dump to matrix market format

            const lduAddressing& lduAddr = matrix_.mesh().lduAddr();
            const labelUList& uPtr = lduAddr.upperAddr();
            const labelUList& lPtr = lduAddr.lowerAddr();
            const labelUList& ownerStart = lduAddr.ownerStartAddr();
            const labelUList& losort = lduAddr.losortAddr();
            const labelUList& losortStart = lduAddr.losortStartAddr();
            const scalarField& diag = matrix_.diag();
            const scalarField& upperPtr = matrix_.upper();
            const scalarField& lowerPtr =
            (
                matrix_.hasLower()
              ? matrix_.lower()
              : matrix_.upper()
            );

            OFstream str("matrixmarket.txt");
            Info<< "Dumping matrix to " << str.name() << endl;

            str << "%%MatrixMarket matrix coordinate real general" << nl
                << lduAddr.size() << ' ' << lduAddr.size() << ' '
                << diag.size()+uPtr.size()+lPtr.size() << nl;

            forAll(diag, celli)
            {
                // Lower
                for
                (
                    label i = losortStart[celli];
                    i < losortStart[celli+1];
                    i++
                )
                {
                    label facei = losort[i];
                    str << celli+1 << ' ' << lPtr[facei]+1 << ' '
                        << lowerPtr[facei] << nl;
                }

                // Diag
                str << celli+1 << ' ' << celli+1 << ' '
                    << diag[celli] << nl;

                // Upper
                for
                (
                    label facei = ownerStart[celli];
                    facei < ownerStart[celli+1];
                    facei++
                )
                {
                    str << celli+1 << ' ' << uPtr[facei]+1 << ' '
                        << upperPtr[facei] << nl;
                }
            }
        }



        if (debug & 2)
        {
            const lduAddressing& lduAddr = matrix_.mesh().lduAddr();
            const labelUList& uPtr = lduAddr.upperAddr();
            const labelUList& lPtr = lduAddr.lowerAddr();
            const scalarField& diag = matrix_.diag();
            const scalarField& upperPtr = matrix_.upper();
            const scalarField& lowerPtr =
            (
                matrix_.hasLower()
              ? matrix_.lower()
              : matrix_.upper()
            );

            Pout<< "Face-based matrix:" << endl;
            forAll(source, celli)
            {
                Pout<< "cell:" << celli
                    << " diag:" << diag[celli]
                    << " source:" << source[celli] << endl;
            }
            forAll(uPtr, facei)
            {
                Pout<< "face:" << facei
                    << " u:" << uPtr[facei] << " coeff:" << upperPtr[facei]
                    << " l:" << lPtr[facei] << " coeff:" << lowerPtr[facei]
                    << endl;
            }
            forAll(interfaceBouCoeffs_, patchi)
            {
                if (interfaceBouCoeffs_.set(patchi))
                {
                    const scalarField& iCoeffs = interfaceIntCoeffs_[patchi];
                    const scalarField& bCoeffs = interfaceBouCoeffs_[patchi];
                    const labelUList& fc = lduAddr.patchAddr(patchi);
                    Pout<< "Patch:" << patchi << endl;
                    forAll(fc, i)
                    {
                        Pout<< "    cell:" << fc[i]
                            << " iCoeff:" << iCoeffs[i]
                            << " bCoeff:" << bCoeffs[i]
                            << endl;
                    }
                }
            }


            Pout<< "CRS matrix:" << endl;
            for (unsigned int celli = 0; celli < rhs.size(); celli++)
            {
                Pout<< "cell:" << celli << " rhs:" << rhs[celli] << endl;
                for (label facei = ptr_[celli]; facei < ptr_[celli+1]; facei++)
                {
                    Pout<< "    cell:" << col[facei]
                        << " coeff:" << val[facei] << endl;
                }
            }
        }


        std::vector<double> solution(rhs.size());
        forAll(psi, celli)
        {
            solution[celli] = psi[celli];
        }


        typedef amgcl::backend::builtin<double> Backend;


        ptrdiff_t n(globalNumbering_.localSize());

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
                amgcl::solver::skyline_lu<double>
            > Solver;

            Solver::params prm;
            prm.isolver.tol = tol;
            prm.num_def_vec = 1;
            prm.def_vec = constant_deflation();
            Backend::params bprm;

            amgcl::mpi::communicator world(MPI_COMM_WORLD);

            ptrdiff_t chunk(n);

//            Solver solve(world, std::tie(chunk, ptr_, col, val), prm, bprm);
//
//             if (debug)
//             {
//                 Pout<< "Constructed AMG solver in "
//                     << solverTime.cpuTimeIncrement() << endl;
//             }
// 
//             std::tie(iters, resid) = solve(rhs, solution);
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
            prm.solver.tol = tol;
            Backend::params bprm;

            Solver solve(std::tie(n, ptr_, col, val), prm, bprm);

            if (debug)
            {
                Pout<< "Constructed AMG solver in "
                    << solverTime.cpuTimeIncrement() << endl;
            }

            std::tie(iters, resid) = solve(rhs, solution);
        }

        if (debug)
        {
            Pout<< "Solved with rhs in " << solverTime.cpuTimeIncrement()
                << endl;
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

        if (false)
        {
            double norm_rhs = sqrt(amgcl::backend::inner_product(rhs, rhs));
            amgcl::backend::spmv
            (
                -1,
                std::tie(n, ptr_, col, val),
                solution,
                1,
                rhs
            );
            double realResid =
                sqrt(amgcl::backend::inner_product(rhs, rhs))
              / norm_rhs;
            DebugVar(realResid);
            DebugVar(resid);
            DebugVar(solverPerf.finalResidual());
        }
    }

    return solverPerf;
}


// ************************************************************************* //
