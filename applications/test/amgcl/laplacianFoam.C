/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "amgcl/amg.hpp"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/adapter/crs_tuple.hpp"
#include "amgcl/coarsening/ruge_stuben.hpp"
#include "amgcl/coarsening/pointwise_aggregates.hpp"
#include "amgcl/coarsening/aggregation.hpp"
#include "amgcl/coarsening/smoothed_aggregation.hpp"
#include "amgcl/relaxation/damped_jacobi.hpp"
#include "amgcl/relaxation/spai0.hpp"
#include "amgcl/solver/gmres.hpp"

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );

            fvOptions.constrain(TEqn);

            //TEqn.solve();

            const labelUList& uPtr = TEqn.lduAddr().upperAddr();
            const labelUList& lPtr = TEqn.lduAddr().lowerAddr();
            //const labelUList& ownStartPtr = TEqn.lduAddr().ownerStartAddr();
            const label n = TEqn.lduAddr().size();

            const scalarList& diag = TEqn.diag();
            const scalarList& upperPtr = TEqn.upper();
            // Note: force symmetry for ease of conversion
            const scalarList& lowerPtr = TEqn.lower();

            // Count number of neighbours
            labelList nNbrs(n, 0);
            forAll(uPtr, facei)
            {
                nNbrs[uPtr[facei]]++;
                nNbrs[lPtr[facei]]++;
            }

            std::vector<double> val(n+lPtr.size()+uPtr.size());   // Values of nonzero entries.
            std::vector<int> col(val.size());       // Column numbers of nonzero entries.
            std::vector<int> ptr(n+1);              // Points to the start of each row in the above arrays.

            ptr[0] = 0;
            for (label celli = 0; celli < n; celli++)
            {
                ptr[celli+1] = ptr[celli]+1+nNbrs[celli];
            }

            std::vector<int> offsets(ptr);
            // Store connections to lower numbered cells
            forAll(uPtr, facei)
            {
                int& index = offsets[uPtr[facei]];
                col[index] = lPtr[facei];
                val[index++] = lowerPtr[facei];
            }
            // Store diagonal
            for (label celli = 0; celli < n; celli++)
            {
                int& index = offsets[celli];
                col[index] = celli;
                val[index++] = diag[celli];
            }
            // Store connections to higher numbered cells
            forAll(lPtr, facei)
            {
                int& index = offsets[lPtr[facei]];
                col[index] = uPtr[facei];
                val[index++] = upperPtr[facei];
            }

            std::vector<double> rhs(n);    // Right-hand side of the system of equations.
            forAll(TEqn.source(), celli)
            {
                rhs[celli] = TEqn.source()[celli];
            }
            lduInterfaceFieldPtrsList interfaces =
                TEqn.psi().boundaryField().scalarInterfaces();
            forAll(TEqn.internalCoeffs(), patchi)
            {
                const scalarField& iCoeffs = TEqn.internalCoeffs()[patchi];
                const scalarField& bCoeffs = TEqn.boundaryCoeffs()[patchi];

                const labelUList& fc = TEqn.lduAddr().patchAddr(patchi);
                forAll(fc, i)
                {
                    rhs[fc[i]] -= iCoeffs[i]*bCoeffs[i];
                }
            }


            typedef amgcl::backend::builtin<double> Backend;

            // Define the AMG type:
            typedef amgcl::amg<
                Backend,
                amgcl::coarsening::smoothed_aggregation,
                amgcl::relaxation::spai0
                > AMG;

            // Construct the AMG hierarchy.
            // Note that this step only depends on the matrix. Hence, the constructed
            // instance may be reused for several right-hand sides.
            // The matrix is specified as a tuple of sizes and ranges.
            AMG::params prm;
            prm.coarse_enough = 10;
            AMG amg( boost::tie(n, ptr, col, val), prm );
    //        AMG amg( boost::tie(n, ptr, col, val) );

            // Output some information about the constructed hierarchy:
            std::cout << amg << std::endl;

            // Use GMRES as an iterative solver:
            typedef amgcl::solver::gmres<Backend> Solver;

            // Construct the iterative solver. It needs size of the system to
            // preallocate the required temporary structures:
            Solver solve(n);

            // Solve the system. Returns number of iterations made and the achieved residual.
            label iters;
            scalar resid;
            std::vector<double>  solution(n);
            boost::tie(iters, resid) = solve(amg, rhs, solution);
            Info << iters << "   " << resid << endl;
            //std::cout<< "solution:" << solution << std::endl;

            fvOptions.correct(T);
        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
