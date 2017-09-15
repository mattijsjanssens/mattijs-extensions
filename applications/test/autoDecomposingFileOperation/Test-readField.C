/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Test-readField

Description
    Read volScalarField

\*---------------------------------------------------------------------------*/

#include "GeometricField.H"
#include "argList.H"
#include "uVolFields.H"
#include "unallocatedFvBoundaryMesh.H"
#include "unallocatedFvMesh.H"
//#include "unallocatedFvPatchField.H"
#include "unallocatedGenericFvPatchField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!Pstream::parRun())
    {
        FatalErrorInFunction << "Not running in parallel" << exit(FatalError);
    }

    // Parent database
    Time baseRunTime
    (
        runTime.controlDict(),
        runTime.rootPath(),
        runTime.globalCaseName(),
        runTime.system(),
        runTime.constant(),
        false                   // enableFunctionObjects
    );
    baseRunTime.setTime(runTime);


    // Extract parent mesh information
    label baseNCells = -1;
    labelList basePatchSizes;
    wordList basePatchNames;
    // TBD: basePatchGroups
    {
        // For now read in base mesh. TBD.

        polyMesh baseMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                baseRunTime.timeName(),
                baseRunTime,
                IOobject::MUST_READ
            )
        );

        baseNCells = baseMesh.nCells();

        const polyBoundaryMesh& pbm = baseMesh.boundaryMesh();
        basePatchNames = pbm.names();
        basePatchSizes.setSize(pbm.size());
        forAll(pbm, patchi)
        {
            basePatchSizes[patchi] = pbm[patchi].size();
        }
    }

    const fvBoundaryMesh& fp = mesh.boundary();

    // Construct an unallocated boundary. (note: the references to the
    // mesh-boundary are not used; names and sizes are from the base mesh)
    unallocatedFvBoundaryMesh baseBoundary;
    baseBoundary.setSize(basePatchNames.size());
    forAll(basePatchNames, patchi)
    {
        baseBoundary.set
        (
            patchi,
            new unallocatedGenericFvPatch
            (
                fp[patchi].patch(),                 // unused reference
                basePatchNames[patchi],
                basePatchSizes[patchi],
                fp                                  // unused reference
            )
        );
    }

    // Construct the mesh
    unallocatedFvMesh baseMesh
    (
        baseRunTime,            // database
        baseNCells,             // nCells
        baseBoundary,           // boundary
        mesh.globalData()       // used for patchSchedule() only
    );


    uVolScalarField p
    (
        IOobject
        (
            "p",
            baseRunTime.timeName(),
            baseMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        baseMesh
    );
    //uVolScalarField p
    //(
    //    IOobject
    //    (
    //        "p",
    //        baseRunTime.timeName(),
    //        baseMesh.thisDb(),
    //        IOobject::NO_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    baseMesh,
    //    dimensionedScalar("zero", dimless, Zero),
    //    unallocatedGenericFvPatchField<scalar>::typeName
    //);
    DebugVar(p);
    p.rename("my_p");
    {
        const bool master = Pstream::master();
        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;
        if (master)
        {
            p.write();
        }
        Pstream::parRun() = oldParRun;
    }

//    volScalarField p
//    (
//        IOobject
//        (
//            "p",
//            runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );
//    volVectorField U
//    (
//        IOobject
//        (
//            "U",
//            runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );

    return 0;
}


// ************************************************************************* //
