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

            if (i > 0)
            {
                const volVectorField& fineU = UFluid[i-1];
                volVectorField fineRes("res", fineU-fineU.prevIter());
                DebugVar(fineRes);

                const meshToMesh& mapper = fluidMappers[i-1];
                tmp<volVectorField> tcoarseRes(mapper.mapSrcToTgt(fineRes));
                DebugVar(tcoarseRes());
                U += tcoarseRes;
            }

            #include "UEqn.H"
            #include "pEqn.H"
            turb.correct();
        }


        Info<< "\nCoarse to fine solution" << endl;
        for (label i = fluidRegions.size()-2; i >= 0; --i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"

            U.storePrevIter();

            {
                const volVectorField& coarseU = UFluid[i+1];
                volVectorField coarseRes("res", coarseU-coarseU.prevIter());
                DebugVar(coarseRes);

                const meshToMesh& mapper = fluidMappers[i];
                tmp<volVectorField> tfineRes(mapper.mapTgtToSrc(coarseRes));
                DebugVar(tfineRes());
                U += tfineRes;
            }

            #include "UEqn.H"
            #include "pEqn.H"
            turb.correct();
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
