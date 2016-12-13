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

DebugVar(pRelax);
DebugVar(URelax);
DebugVar(phiRelax);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Fine to coarse

        Info<< "\nFine to coarse solution" << endl;
        forAll(fluidRegions, i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"

            U.storePrevIter();
            p.storePrevIter();

            if (runTime.outputTime())
            {
                volVectorField("UStart", U).write();
                volScalarField("pStart", p).write();
            }

            // Update U,phi
            if (i > 0)
            {
                const meshToMesh0& mapper = fineToCoarseMappers[i-1];

                const volVectorField& fineU = UFluid[i-1];
                volVectorField fineRes
                (
                    "UrestrictFineRes",
                    fineU-fineU.prevIter()
                );
                if (runTime.outputTime())
                {
                    fineRes.write();
                }

                tmp<volVectorField> tcoarseRes
                (
                    mapper.interpolate
                    (
                        fineRes,
                        meshToMesh0::INTERPOLATE,
                        eqOp<vector>()
                    )
                );

                //DebugVar(tcoarseRes());
                if (runTime.outputTime())
                {
                    tcoarseRes.ref().rename("UrestrictCoarseRes");
                    tcoarseRes().write();
                }

                phi += phiRelax*fvc::flux(tcoarseRes());
                U += URelax*tcoarseRes;

                if (runTime.outputTime())
                {
                    volVectorField("UPredict", U).write();
                }
            }
            // Update p
            if (i > 0)
            {
                const volScalarField& fineP = pFluid[i-1];
                volScalarField fineRes
                (
                    "prestrictFineRes",
                    fineP-fineP.prevIter()
                );

                const meshToMesh0& mapper = fineToCoarseMappers[i-1];
                tmp<volScalarField> tcoarseRes
                (
                    mapper.interpolate
                    (
                        fineRes,
                        meshToMesh0::INTERPOLATE,
                        eqOp<scalar>()
                    )
                );

                p += pRelax*tcoarseRes;

                if (runTime.outputTime())
                {
                    fineRes.write();
                    volScalarField("pPredict", p).write();
                }
            }

            #include "UEqn.H"
            #include "pEqn.H"
            turb.correct();

            if (runTime.outputTime())
            {
                volVectorField("UDown", U).write();
                volScalarField("pDown", p).write();
            }
        }


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
                const volVectorField& coarseU = UFluid[i+1];
                volVectorField coarseRes
                (
                    "prolongCoarseRes",
                    coarseU-coarseU.prevIter()
                );

                if (runTime.outputTime())
                {
                    coarseRes.write();
                }

                const meshToMesh0& mapper = coarseToFineMappers[i];
                tmp<volVectorField> tfineRes
                (
                    mapper.interpolate
                    (
                        coarseRes,
                        meshToMesh0::INTERPOLATE,
                        eqOp<vector>()
                    )
                );

                if (runTime.outputTime())
                {
                    tfineRes.ref().rename("prolongFineRes");
                    tfineRes().write();
                }

                phi += phiRelax*fvc::flux(tfineRes());
                U += URelax*tfineRes;

                if (runTime.outputTime())
                {
                    volVectorField("U_updated_from_coarse", U).write();
                }
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

                if (runTime.outputTime())
                {
                    volScalarField("p_updated_from_coarse", p).write();
                }
            }

            #include "UEqn.H"
            #include "pEqn.H"
            turb.correct();

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
