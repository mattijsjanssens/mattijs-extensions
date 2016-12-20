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
//#include "meshToMesh.H"
#include "meshToMesh0.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void interpolate
(
    Field<Type>& dest,
    const meshToMesh0& mapper,
    const Field<Type>& source
)
{
    // Helper function to interpolate from Field to Field

    GeometricField<Type, fvPatchField, volMesh> res
    (
        IOobject
        (
            "res",
            mapper.fromMesh().time().timeName(),
            mapper.fromMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mapper.fromMesh(),
        dimensioned<Type>
        (
            "zero",
            dimless,
            pTraits<Type>::zero
        )
    );
    res.primitiveFieldRef() = source;

    mapper.interpolateInternalField
    (
        dest,
        res,
        meshToMesh0::INTERPOLATE,
        eqOp<Type>()
    );
}


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

            // Update the defect
            if (i > 0)
            {
                #include "assembleUEqn.H"
                tmp<fvVectorMatrix> tUpEqn(UEqn == -fvc::grad(p));
                UDefect = URes;
                UDefect.Field<vector>::operator-=(tUpEqn().residual());
                //pDefect = pRes;
                //pDefect -= fvc::div(phi)();
            }


            if (i == fluidRegions.size()-1)
            {
                p.storePrevIter();
                U.storePrevIter();
            }


            // Restrict to coarse mesh

            Pout<< "pRes:" << gAverage(pRes)
                << " pDefect:" << gAverage(pDefect)
                << "  URes:" << gAverage(URes)
                << " UDefect:" << gAverage(UDefect)
                << endl;


            label nSmooth = 0;  //(i == 0 ? 1: 4);
            for (label smoothi = 0; smoothi < nSmooth; smoothi++)
            {
                #include "UEqn.H"
                #include "pEqn.H"

                transport.correct();
                turb.correct();
            }



            // Re-evaluate residual and map to coarser level
            if (i < fluidRegions.size()-1)
            {
                #include "assembleUEqn.H"

                tmp<fvVectorMatrix> tUpEqn(UEqn == -fvc::grad(p));
                interpolate
                (
                    UFluidRes[i+1],
                    fineToCoarseMappers[i],
                    tUpEqn().residual()()
                );

                interpolate
                (
                    pFluidRes[i+1],
                    fineToCoarseMappers[i],
                    fvc::div(phi)()
                );
            }
        }


        // Prolongation - coarse to fine

        Info<< "\nCoarse to fine solution" << endl;
        for (label i = fluidRegions.size()-2; i >= 1; --i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"

            U.storePrevIter();
            p.storePrevIter();


            // Update fields
            {
                const meshToMesh0& mapper = coarseToFineMappers[i];

                const volVectorField& coarseU = UFluid[i+1];
                tmp<volVectorField> tfineDeltaU
                (
                    mapper.interpolate
                    (
                        coarseU-coarseU.prevIter(),
                        meshToMesh0::INTERPOLATE,
                        eqOp<vector>()
                    )
                );

                phi += phiRelax*fvc::flux(tfineDeltaU());
                U += URelax*tfineDeltaU;


                const volScalarField& coarseP = pFluid[i+1];
                tmp<volScalarField> tfineDeltaP
                (
                    mapper.interpolate
                    (
                        coarseP-coarseP.prevIter(),
                        meshToMesh0::INTERPOLATE,
                        eqOp<scalar>()
                    )
                );

                p += pRelax*tfineDeltaP;
            }


            label nSmooth = 1;
            for (label smoothi = 0; smoothi < nSmooth; smoothi++)
            {
                #include "UEqn.H"
                #include "pEqn.H"
                transport.correct();
                turb.correct();
            }

            if (i == 1)
            {
                // Map to finest level. Gets solved in next iteration.

                const meshToMesh0& mapper = coarseToFineMappers[0];

                tmp<volVectorField> tfineDeltaU
                (
                    mapper.interpolate
                    (
                        U-U.prevIter(),
                        meshToMesh0::INTERPOLATE,
                        eqOp<vector>()
                    )
                );
                phiFluid[0] += phiRelax*fvc::flux(tfineDeltaU());
                UFluid[0] += URelax*tfineDeltaU;


                tmp<volScalarField> tfineDeltaP
                (
                    mapper.interpolate
                    (
                        p-p.prevIter(),
                        meshToMesh0::INTERPOLATE,
                        eqOp<scalar>()
                    )
                );
                pFluid[0] += pRelax*tfineDeltaP;
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
