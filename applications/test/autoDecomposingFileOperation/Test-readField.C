/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
#include "parFvFieldReconstructor.H"
#include "polyBoundaryMeshEntries.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<mapDistribute> calcReconstructMap
(
    const labelList& localToGlobal,
    const bool localHasFlip
)
{
    // localToGlobal: which elements we want from the master and in what order

    label globalSize;
    if (localHasFlip)
    {
        globalSize = gMax(mag(localToGlobal)-1)+1;
    }
    else
    {
        globalSize = gMax(localToGlobal)+1;
    }

    // On master: receive all the undecomposed elements
    labelListList subMap(Pstream::nProcs());
    subMap[Pstream::myProcNo()] = localToGlobal;
    Pstream::gatherList(subMap);
    if (!Pstream::master())
    {
        // No need to receive anything on slaves
        subMap = labelList(0);
    }

    // Everywhere: send subset to master
    labelListList constructMap(Pstream::nProcs());
    constructMap[0] = identity(localToGlobal.size());

    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            globalSize,             // max index of receive data
            constructMap.xfer(),    // what to send
            subMap.xfer(),          // what to receive
            false,                  // no sign flag in sending
            localHasFlip            // have sign in receiving
        )
    );
}


autoPtr<mapDistributePolyMesh> readProcMap
(
    const IOobject& io,
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    const labelList& oldPatchStarts,
    const labelList& oldPatchNMeshPoints
)
{
    // Read local cell allocation
    labelIOList cellMap(IOobject(io, "cellProcAddressing"));
    autoPtr<mapDistribute> cellMapPtr(calcReconstructMap(cellMap, false));

    labelIOList faceMap(IOobject(io, "faceProcAddressing"));
    autoPtr<mapDistribute> faceMapPtr(calcReconstructMap(faceMap, true));

    labelIOList pointMap(IOobject(io, "pointProcAddressing"));
    autoPtr<mapDistribute> pointMapPtr(calcReconstructMap(pointMap, false));

    labelIOList boundaryMap(IOobject(io, "boundaryProcAddressing"));
    autoPtr<mapDistribute> bMapPtr(calcReconstructMap(boundaryMap, false));
DebugVar(bMapPtr());

    return autoPtr<mapDistributePolyMesh>
    (
        new mapDistributePolyMesh
        (
            nOldPoints,
            nOldFaces,
            nOldCells,
            xferCopy(oldPatchStarts),
            xferCopy(oldPatchNMeshPoints),

            // how to subset pieces of mesh to send across
            pointMapPtr().xfer(),
            faceMapPtr().xfer(),
            cellMapPtr().xfer(),
            bMapPtr().xfer()
        )
    );
}


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


    // Read procAddressing files. Deduct base mesh sizes.
    autoPtr<mapDistributePolyMesh> distMap
    (
        readProcMap
        (
            IOobject
            (
                "cellProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            0,              //nOldPoints
            0,              //nOldFaces
            0,              //nOldCells
            labelList(0),   //oldPatchStarts
            labelList(0)    //oldPatchNMeshPoints
        )
    );


    // Extract parent mesh information
    labelList basePatchSizes;
    labelList basePatchStarts;
    wordList basePatchNames;
    {
        polyBoundaryMeshEntries patchEntries
        (
            IOobject
            (
                "boundary",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                baseRunTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        basePatchNames.setSize(patchEntries.size());
        basePatchSizes.setSize(patchEntries.size());
        basePatchStarts.setSize(patchEntries.size());
        forAll(patchEntries, patchi)
        {
            basePatchNames[patchi] = patchEntries[patchi].keyword();
            const dictionary& d = patchEntries[patchi].dict();
            basePatchSizes[patchi] = readLabel(d.lookup("nFaces"));
            basePatchStarts[patchi] = readLabel(d.lookup("startFace"));
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
                basePatchStarts[patchi],
                fp                                  // unused reference
            )
        );
    }

    // Construct the (unallocated) mesh
    unallocatedFvMesh baseMesh
    (
        baseRunTime,                            // database
        distMap().cellMap().constructSize(),    // nCells
        baseBoundary,                           // boundary
        mesh.globalData()                       // used for patchSchedule() only
    );

    // // Read field on baseMesh
    // uVolScalarField baseFld
    // (
    //     IOobject
    //     (
    //         "p",
    //         baseRunTime.timeName(),
    //         baseMesh.thisDb(),
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     baseMesh
    // );
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

    parFvFieldReconstructor reconstructor
    (
        baseMesh,
        mesh,
        distMap(),
        false           // isWriteProc
    );

    // Load local field
    volScalarField p
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh
    );

    // Map local field onto baseMesh
    tmp<uVolScalarField> tfld
    (
        reconstructor.reconstructFvVolumeField(p)
    );


    // Write master field to parent
    DebugVar(p);
    //DebugVar(tfld());

    tfld.ref().rename("my_p");
    {
        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;
        if (Pstream::master())
        {
            tfld().write();
        }
        Pstream::parRun() = oldParRun;
    }

    return 0;
}


// ************************************************************************* //
