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
#include "amgcl/make_solver.hpp"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/adapter/crs_tuple.hpp"
#include "amgcl/coarsening/ruge_stuben.hpp"
#include "amgcl/coarsening/pointwise_aggregates.hpp"
#include "amgcl/coarsening/aggregation.hpp"
#include "amgcl/coarsening/smoothed_aggregation.hpp"
#include "amgcl/relaxation/damped_jacobi.hpp"
#include "amgcl/relaxation/spai0.hpp"

#include "amgcl/solver/bicgstab.hpp"
#include "amgcl/solver/gmres.hpp"

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void printMatrix(fvMatrix<Type>& m)
{
    Pout<< "m.psi:" << m.psi() << endl;
    Pout<< "m.diag:" << m.diag() << endl;
    Pout<< "m.source:" << m.source() << endl;
    if (m.hasLower())
    {
        Pout<< "m.lower:" << m.lower() << endl;
    }
    if (m.hasUpper())
    {
        Pout<< "m.upper:" << m.upper() << endl;
    }
    const fvMesh& mesh = m.psi().mesh();
    forAll(m.internalCoeffs(), patchI)
    {
        if (m.internalCoeffs().set(patchI))
        {
            Pout<< "patch:" << mesh.boundaryMesh()[patchI].name() << nl
                << "    internal:" << m.internalCoeffs()[patchI] << nl
                << "    boundary:" << m.boundaryCoeffs()[patchI] << endl; 
        }
    }
    Pout<< endl;
}


void printMatrix
(
    const std::vector<int>& ptr,
    const std::vector<int>& col,
    const std::vector<double>& val,
    const std::vector<double>& rhs
)
{
    for (unsigned int celli = 0; celli < rhs.size(); celli++)
    {
        Pout<< "cell:" << celli << endl;
        DebugVar(rhs[celli]);

        label start = ptr[celli];
        label end = ptr[celli+1];
        for (label facei = start; facei < end; facei++)
        {
            Pout<< "    " << facei << "\tval:" << val[facei]
                << "\tto cell:" << col[facei]
                << endl;
        }
        Pout<< endl;
    }
}


typedef amgcl::backend::builtin<double> Backend;

// // Define the AMG type:
// typedef amgcl::amg
// <
//     Backend,
//     amgcl::coarsening::smoothed_aggregation,
//     amgcl::relaxation::spai0
// > AMG;

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


// void convertWithAddedBoundaryNodes
// (
//     fvScalarMatrix& TEqn,
//     std::vector<double>& val,
//     std::vector<int>& col,
//     std::vector<int>& ptr,
//     std::vector<double>& rhs
// )
// {
//     // Convert face-based, diagonal+off diagonal notation to general sparse
//     // matrix. Adds fixed value bcs
// 
//     const lduAddressing& lduAddr = TEqn.lduAddr();
//     const labelUList& uPtr = lduAddr.upperAddr();
//     const labelUList& lPtr = lduAddr.lowerAddr();
//     const label n = lduAddr.size();
//     const scalarList& diag = TEqn.diag();
//     const scalarList& upperPtr = TEqn.upper();
//     // Note: force symmetry for ease of conversion
//     const scalarList& lowerPtr = TEqn.lower();
// 
// 
//     // Count number of neighbours
//     labelList nNbrs(n, 0);
//     forAll(uPtr, facei)
//     {
//         nNbrs[uPtr[facei]]++;
//         nNbrs[lPtr[facei]]++;
//     }
// 
//     label nBnd = 0;
//     labelList nBndPerCell(n, 0);
//     forAll(TEqn.internalCoeffs(), patchi)
//     {
//         const scalarField& iCoeffs = TEqn.internalCoeffs()[patchi];
//         const scalarField& bCoeffs = TEqn.boundaryCoeffs()[patchi];
// 
//         const labelUList& fc = lduAddr.patchAddr(patchi);
// 
//         Pout<< "patch:" << patchi
//             << " fc:" << fc
//             << " iCoeffs:" << iCoeffs
//             << " bCoeffs:" << bCoeffs
//             << endl;
// 
//         nBnd += fc.size();
//         forAll(fc, i)
//         {
//             nBndPerCell[fc[i]]++;
//         }
//     }
// 
// Pout<< "nCells      :" << n << endl;
// Pout<< "nBnd        :" << nBnd << endl;
// Pout<< "nNbrs       :" << nNbrs << endl;
// Pout<< "nBndPerCell :" << nBndPerCell << endl;
// 
//     rhs.resize(n+nBnd, 0.0);   // Right-hand side
// 
//     ptr.resize(n+nBnd+1);   // Points to the start of each row
//     ptr[0] = 0;
//     for (label celli = 0; celli < n; celli++)
//     {
//         ptr[celli+1] = ptr[celli]+1+nNbrs[celli]+nBndPerCell[celli];
//     }
//     for (label celli = n; celli < n+nBnd; celli++)
//     {
//         ptr[celli+1] = ptr[celli]+1;
//     }
// 
//     col.resize(ptr[ptr.size()-1]);
//     val.resize(ptr[ptr.size()-1]);
// 
// 
//     std::vector<int> offsets(ptr);
//     // Store connections to lower numbered cells
//     forAll(uPtr, facei)
//     {
//         int& index = offsets[uPtr[facei]];
//         col[index] = lPtr[facei];
//         val[index++] = lowerPtr[facei];
//     }
//     // Store diagonal
//     for (label celli = 0; celli < n; celli++)
//     {
//         int& index = offsets[celli];
//         col[index] = celli;
//         val[index++] = diag[celli];
//     }
//     // Store connections to higher numbered cells
//     forAll(lPtr, facei)
//     {
//         int& index = offsets[lPtr[facei]];
//         col[index] = uPtr[facei];
//         val[index++] = upperPtr[facei];
//     }
//     // Store connections to boundary values
//     label bndCelli = n;
//     forAll(TEqn.internalCoeffs(), patchi)
//     {
//         const labelUList& fc = TEqn.lduAddr().patchAddr(patchi);
//         forAll(fc, i)
//         {
//             int& index = offsets[fc[i]];
//             col[index] = bndCelli;
//             val[index++] = -TEqn.internalCoeffs()[patchi][i];
// 
//             // Do boundary cell itself
//             {
//                 index = ptr[bndCelli];
//                 col[index] = bndCelli;
//                 val[index] = 1.0;   //TEqn.boundaryCoeffs()[patchi][i];
//                 rhs[bndCelli] = TEqn.boundaryCoeffs()[patchi][i];
//             }
//             bndCelli++;
//         }
//     }    
// 
//     // Source for cells (not boundary cells)
//     forAll(TEqn.source(), celli)
//     {
//         rhs[celli] = TEqn.source()[celli];
//     }
// }


void convert
(
    fvScalarMatrix& TEqn,
    std::vector<double>& val,
    std::vector<int>& col,
    std::vector<int>& ptr,
    std::vector<double>& rhs
)
{
    // Convert face-based, diagonal+off diagonal notation to general sparse
    // matrix

    const lduAddressing& lduAddr = TEqn.lduAddr();
    const labelUList& uPtr = lduAddr.upperAddr();
    const labelUList& lPtr = lduAddr.lowerAddr();
    const label n = lduAddr.size();
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

    // Values of nonzero entries
    val.resize(n+lPtr.size()+uPtr.size());
    // Column numbers of nonzero entries
    col.resize(val.size());
    // Points to the start of each row in the above arrays
    ptr.resize(n+1);

    ptr[0] = 0;
    for (label celli = 0; celli < n; celli++)
    {
        // Space for diagonal and cell neighbours
        ptr[celli+1] = ptr[celli]+1+nNbrs[celli];
    }

    std::vector<int> offsets(ptr);
    // Store diagonal
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

    rhs.resize(n); // Right-hand side of the system of equations
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

        //Pout<< "patch:" << patchi << " fc:" << fc
        //    << " iCoeffs:" << iCoeffs << " bCoeffs:" << bCoeffs << endl;

        forAll(fc, i)
        {
            label celli = fc[i];

            label index = ptr[celli];

            if (col[index] != celli)
            {
                FatalErrorInFunction << "Diagonal not first" << exit(FatalError);
            }

            val[index] += iCoeffs[i];
            rhs[celli] += bCoeffs[i];
        }
    }
}


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

            Info<< "After fv discretisation = "
                << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            // Force symmetry for ease of conversion
            (void)TEqn.lower();

            //printMatrix(TEqn);

            std::vector<double> val;
            std::vector<int> col;
            std::vector<int> ptr;
            std::vector<double> rhs;
            convert(TEqn, val, col, ptr, rhs);
            const label n = rhs.size();

            Info<< "After conversion = "
                << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            // Construct the AMG hierarchy.
            // Note that this step only depends on the matrix. Hence,
            // the constructed
            // instance may be reused for several right-hand sides.
            // The matrix is specified as a tuple of sizes and ranges.
            //AMG::params prm;
            //prm.coarse_enough = 10;

            //AMG amg( boost::tie(n, ptr, col, val), prm );
            //AMG amg( boost::tie(n, ptr, col, val) );

            // Output some information about the constructed hierarchy:
            //std::cout << amg << std::endl;

            // Use GMRES as an iterative solver:
            //typedef amgcl::solver::gmres<Backend> Solver;

            // Construct the iterative solver. It needs size of the system to
            // preallocate the required temporary structures:
            Solver::params prm;
            prm.solver.tol = 1e-6;
            Solver solve(boost::tie(n, ptr, col, val), prm);

            Info<< "After constructing solver = "
                << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            // Solve the system. Returns number of iterations made and the
            // achieved residual.
            label iters;
            scalar resid;
            std::vector<double>  solution(n, 0.0);
            boost::tie(iters, resid) = solve(rhs, solution);
            Info << iters << "   " << resid << endl;

            Info<< "After AMGCL solving = "
                << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            forAll(solution, celli)
            {
                T[celli] = solution[celli];
            }
            T.correctBoundaryConditions();
            Info<< "AMGCL residual:" << gAverage(TEqn.residual()) << endl;


            //std::cout<< "solution:" << solution << std::endl;

            //for (unsigned int celli = 0; celli < solution.size(); celli++)
            //{
            //    DebugVar(solution[celli]);
            //}

            T = dimensionedScalar("zero", T.dimensions(), 0.0);
            TEqn.solve();

            Info<< "After OpenFOAM solving = "
                << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            Info<< "OpenFOAM residual:" << gAverage(TEqn.residual()) << endl;

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
