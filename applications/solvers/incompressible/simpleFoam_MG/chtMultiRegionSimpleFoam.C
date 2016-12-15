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

Application
    chtMultiRegionSimpleFoam

Description
    Steady-state solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "regionProperties.H"
#include "meshToMesh.H"
#include "meshToMesh0.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
//     #define NO_CONTROL
//     #define CREATE_MESH createMeshesPostProcess.H
//     #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
//    #include "createControl.H"
    simpleControl simple(fluidRegions[0]);


    Info<< "\nStarting time loop\n" << endl;

    const scalar pRelax = 1;
    const scalar URelax = 1;
    const scalar phiRelax = 1.0;    //0.0;
//
//DebugVar(pRelax);
//DebugVar(URelax);
//DebugVar(phiRelax);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Restriction - fine to coarse

        Info<< "\nFine to coarse solution" << endl;
        forAll(fluidRegions, i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"


            // Restrict to coarse mesh source terms
            if (i > 0)
            {
                const meshToMesh0& mapper = fineToCoarseMappers[i-1];

                tmp<volVectorField> tcoarseResU
                (
                    mapper.interpolate
                    (
                        UFluidRes[i-1],
                        meshToMesh0::INTERPOLATE,
                        eqOp<vector>()
                    )
                );
                tmp<volScalarField> tcoarseResP
                (
                    mapper.interpolate
                    (
                        pFluidRes[i-1],
                        meshToMesh0::INTERPOLATE,
                        eqOp<scalar>()
                    )
                );

                label nSmooth = 1;  //(i == 0 ? 1: 4);
                for (label smoothi = 0; smoothi < nSmooth; smoothi++)
                {
                    // UEqn.H
                    MRF.correctBoundaryVelocity(U);

                    tmp<fvVectorMatrix> tUEqn
                    (
                        fvm::div(phi, U)
                      + MRF.DDt(U)
                      + turb.divDevReff(U)
                     ==
                        fvOptions(U) + tcoarseResU
                    );
                    fvVectorMatrix& UEqn = tUEqn.ref();

                    UEqn.relax();

                    fvOptions.constrain(UEqn);

                    if (simple.momentumPredictor())
                    {
                        tmp<fvVectorMatrix> UpEqn(UEqn == -fvc::grad(p));
                        fvVectorMatrix& m = UpEqn.ref();
                        solve(m);

                        URes.primitiveFieldRef() = m.residual();

                        fvOptions.correct(U);
                    }


                    // pEqn.H
                    {
                        volScalarField rAU(1.0/UEqn.A());
                        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
                        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
                        MRF.makeRelative(phiHbyA);
                        adjustPhi(phiHbyA, U, p);

                        tmp<volScalarField> rAtU(rAU);

                        if (simple.consistent())
                        {
                            rAtU = 1.0/(1.0/rAU - UEqn.H1());
                            phiHbyA +=
                                fvc::interpolate(rAtU() - rAU)
                               *fvc::snGrad(p)*mesh.magSf();
                            HbyA -= (rAU - rAtU())*fvc::grad(p);
                        }

                        tUEqn.clear();

                        // Update the pressure BCs to ensure flux consistency
                        constrainPressure(p, U, phiHbyA, rAtU(), MRF);

                        // Non-orthogonal pressure corrector loop
                        while (simple.correctNonOrthogonal())
                        {
                            fvScalarMatrix pEqn
                            (
                                fvm::laplacian(rAtU(), p)
                             ==
                                fvc::div(phiHbyA)
                              + tcoarseResP
                            );

                            pEqn.setReference(pRefCell, pRefValue);

                            pEqn.solve();
                            pRes.primitiveFieldRef() = pEqn.residual();

                            if (simple.finalNonOrthogonalIter())
                            {
                                phi = phiHbyA - pEqn.flux();
                            }
                        }

                        #include "continuityErrs.H"

                        // Explicitly relax pressure for momentum corrector
                        p.relax();

                        // Momentum corrector
                        U = HbyA - rAtU()*fvc::grad(p);
                        U.correctBoundaryConditions();
                        fvOptions.correct(U);
                    }
                }

                turb.correct();
            }
        }


        // Prolongation - fine to coarse

        Info<< "\nCoarse to fine solution" << endl;
        for (label i = fluidRegions.size()-2; i >= 1; --i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"

            U.storePrevIter();
            p.storePrevIter();


            // U, phi
            {
                const meshToMesh0& mapper = coarseToFineMappers[i];

                const volVectorField& coarseU = UFluid[i+1];
                volVectorField coarseRes
                (
                    "prolongCoarseRes",
                    coarseU-coarseU.prevIter()
                );

                tmp<volVectorField> tfineRes
                (
                    mapper.interpolate
                    (
                        coarseRes,
                        meshToMesh0::INTERPOLATE,
                        eqOp<vector>()
                    )
                );

                phi += phiRelax*fvc::flux(tfineRes());
                U += URelax*tfineRes;
            }

            // p
            {
                const volScalarField& coarseP = pFluid[i+1];
                volScalarField coarseRes
                (
                    "prolongCoarseRes",
                    coarseP-coarseP.prevIter()
                );

                const meshToMesh0& mapper = coarseToFineMappers[i];
                tmp<volScalarField> tfineRes
                (
                    mapper.interpolate
                    (
                        coarseRes,
                        meshToMesh0::INTERPOLATE,
                        eqOp<scalar>()
                    )
                );

                p += pRelax*tfineRes;
            }


            label nSmooth = 1;
            for (label smoothi = 0; smoothi < nSmooth; smoothi++)
            {
                #include "UEqn.H"
                #include "pEqn.H"
                turb.correct();
            }

            if (i == 1)
            {
                // Map to finest level. Gets solved in next iteration.

                const meshToMesh0& mapper = coarseToFineMappers[0];

                tmp<volVectorField> tURes
                (
                    mapper.interpolate
                    (
                        U-U.prevIter(),
                        meshToMesh0::INTERPOLATE,
                        eqOp<vector>()
                    )
                );
                phiFluid[0] += phiRelax*fvc::flux(tURes());
                UFluid[0] += URelax*tURes;


                tmp<volScalarField> tpRes
                (
                    mapper.interpolate
                    (
                        p-p.prevIter(),
                        meshToMesh0::INTERPOLATE,
                        eqOp<scalar>()
                    )
                );
                pFluid[0] += pRelax*tpRes;
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
