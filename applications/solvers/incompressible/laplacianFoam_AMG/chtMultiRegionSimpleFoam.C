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
    chtMultiRegionSimpleFoam

Description
    Steady-state solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "rhoThermo.H"
//#include "turbulentFluidThermoModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
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
//     #include "initContinuityErrs.H"


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Fine to coarse

        Info<< "\nFine to coarse solution" << endl;
        forAll(solidRegions, i)
        {
            Info<< "\nSolving for solid region "
                << solidRegions[i].name() << endl;
            #include "setRegionSolidFields.H"

            T.storePrevIter();

            if (i > 0)
            {
                const volScalarField& fineT = Ts[i-1];
                volScalarField fineRes("res", fineT-fineT.prevIter());
                DebugVar(fineRes);

                const meshToMesh& mapper = solidMappers[i-1];
                tmp<volScalarField> tcoarseRes(mapper.mapSrcToTgt(fineRes));
                DebugVar(tcoarseRes());
                T += tcoarseRes;
            }


            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }


        Info<< "\nCoarse to fine solution" << endl;
        for (label i = solidRegions.size()-2; i >= 0; --i)
        {
            Info<< "\nSolving for solid region "
                << solidRegions[i].name() << endl;
            #include "setRegionSolidFields.H"

            T.storePrevIter();

            {
                const volScalarField& coarseT = Ts[i+1];
                volScalarField coarseRes("res", coarseT-coarseT.prevIter());
                DebugVar(coarseRes);

                const meshToMesh& mapper = solidMappers[i];
                tmp<volScalarField> tfineRes(mapper.mapTgtToSrc(coarseRes));
                DebugVar(tfineRes());
                T += tfineRes;
            }

            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
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
