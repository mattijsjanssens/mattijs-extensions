/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "unallocatedFvMeshTools.H"
#include "polyBoundaryMeshEntries.H"
//#include "unallocatedFvBoundaryMesh.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistribute>
Foam::unallocatedFvMeshTools::calcReconstructMap
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


Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::unallocatedFvMeshTools::readReconstructMap(const IOobject& io)
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

    labelList oldPatchStarts;
    labelList oldPatchNMeshPoints;
    return autoPtr<mapDistributePolyMesh>
    (
        new mapDistributePolyMesh
        (
            0,              //nOldPoints,
            0,              //nOldFaces,
            0,              //nOldCells,
            oldPatchStarts.xfer(),          //xferCopy(oldPatchStarts),
            oldPatchNMeshPoints.xfer(),     //xferCopy(oldPatchNMeshPoints),

            // how to subset pieces of mesh to send across
            pointMapPtr().xfer(),
            faceMapPtr().xfer(),
            cellMapPtr().xfer(),
            bMapPtr().xfer()
        )
    );
}


//Foam::autoPtr<Foam::unallocatedFvMesh> Foam::unallocatedFvMeshTools::newMesh
//(
//    const IOobject& io,
//    const fvMesh& mesh,
//    const label nCells
//)
//{
//    // Extract parent mesh information
//    labelList basePatchSizes;
//    labelList basePatchStarts;
//    wordList basePatchNames;
//    wordList basePatchTypes;
//    {
//        IOobject meshIO
//        (
//            io.name(),              // name of mesh
//            io.instance(),
//            fvMesh::meshSubDir,     // location of mesh files
//            io.db(),
//            IOobject::MUST_READ
//        );
//        polyBoundaryMeshEntries patchEntries(IOobject(meshIO, "boundary"));
//
//        basePatchNames.setSize(patchEntries.size());
//        basePatchTypes.setSize(patchEntries.size());
//        basePatchSizes.setSize(patchEntries.size());
//        basePatchStarts.setSize(patchEntries.size());
//        forAll(patchEntries, patchi)
//        {
//            basePatchNames[patchi] = patchEntries[patchi].keyword();
//            const dictionary& d = patchEntries[patchi].dict();
//            basePatchTypes[patchi] = d.lookup("type");
//            basePatchSizes[patchi] = readLabel(d.lookup("nFaces"));
//            basePatchStarts[patchi] = readLabel(d.lookup("startFace"));
//        }
//    }
//
//
//    const fvBoundaryMesh& fp = mesh.boundary();
//
//    // Construct an unallocated boundary. (note: the references to the
//    // mesh-boundary are not used; names and sizes are from the base mesh)
//    unallocatedFvBoundaryMesh baseBoundary;
//    baseBoundary.setSize(basePatchNames.size());
//    forAll(basePatchNames, patchi)
//    {
//        if (patchi < fp.size())
//        {
//            baseBoundary.set
//            (
//                patchi,
//                new unallocatedGenericFvPatch
//                (
//                    fp[patchi].patch(),                 // unused reference
//                    basePatchNames[patchi],
//                    basePatchTypes[patchi],
//                    basePatchSizes[patchi],
//                    basePatchStarts[patchi],
//                    fp                                  // unused reference
//                )
//            );
//        }
//        else
//        {
//            // Patch non-existing in the boundaryMesh
//            // Make up something
//            baseBoundary.set
//            (
//                patchi,
//                new unallocatedGenericFvPatch
//                (
//                    fp.last().patch(),                 // unused reference
//                    basePatchNames[patchi],
//                    basePatchSizes[patchi],
//                    basePatchStarts[patchi],
//                    fp                                  // unused reference
//                )
//            );
//        }
//    }
//
//    // Construct the (unallocated) mesh
//    return autoPtr<unallocatedFvMesh>
//    (
//        new unallocatedFvMesh
//        (
//            io.db(),                    // database
//            nCells,                     // nCells
//            baseBoundary,               // boundary
//            mesh.globalData()           // used for patchSchedule() only
//        )
//    );
//}


Foam::autoPtr<Foam::unallocatedFvMesh> Foam::unallocatedFvMeshTools::newMesh
(
    const IOobject& io,
    const label nCells
)
{
    // Extract parent mesh information
    labelList basePatchSizes;
    labelList basePatchStarts;
    wordList basePatchNames;
    wordList basePatchTypes;
    List<wordList> basePatchGroups;
    {
        fileName meshDir
        (
            io.name() == fvMesh::defaultRegion
          ? fileName(fvMesh::meshSubDir)
          : fileName(io.name())/fvMesh::meshSubDir
        );

        IOobject meshIO
        (
            io.name(),              // name of mesh
            io.time().findInstance(meshDir, "faces"), //io.instance(),
            fvMesh::meshSubDir,     // location of mesh files
            io.db(),
            IOobject::MUST_READ
        );

        polyBoundaryMeshEntries patchEntries(IOobject(meshIO, "boundary"));

        basePatchNames.setSize(patchEntries.size());
        basePatchTypes.setSize(patchEntries.size());
        basePatchSizes.setSize(patchEntries.size());
        basePatchStarts.setSize(patchEntries.size());
        basePatchGroups.setSize(patchEntries.size());
        forAll(patchEntries, patchi)
        {
            basePatchNames[patchi] = patchEntries[patchi].keyword();
            const dictionary& d = patchEntries[patchi].dict();
            d.lookup("type") >> basePatchTypes[patchi];
            basePatchSizes[patchi] = readLabel(d.lookup("nFaces"));
            basePatchStarts[patchi] = readLabel(d.lookup("startFace"));
            d.readIfPresent("inGroups", basePatchGroups[patchi]);
        }
    }

    // Construct the (unallocated) mesh
    return autoPtr<unallocatedFvMesh>
    (
        new unallocatedFvMesh
        (
            io.name(),
            io.db(),                    // database
            nCells,                     // nCells
            basePatchNames,
            basePatchTypes,
            basePatchSizes,
            basePatchStarts,
            basePatchGroups,
            *reinterpret_cast<const globalMeshData*>(0)
        )
    );
}


Foam::autoPtr<Foam::unallocatedFvMesh> Foam::unallocatedFvMeshTools::newMesh
(
    const IOobject& io
)
{
    // Extract mesh information
    label nCells;
    labelList basePatchSizes;
    labelList basePatchStarts;
    wordList basePatchNames;
    wordList basePatchTypes;
    List<wordList> basePatchGroups;
    {
        fileName meshDir
        (
            io.name() == fvMesh::defaultRegion
          ? fileName(fvMesh::meshSubDir)
          : fileName(io.name())/fvMesh::meshSubDir
        );

        IOobject meshIO
        (
            io.name(),              // name of mesh
            io.time().findInstance(meshDir, "faces"), //io.instance(),
            fvMesh::meshSubDir,     // location of mesh files
            io.db(),
            IOobject::MUST_READ
        );

        labelIOList owner(IOobject(meshIO, "owner"));
        if (owner.size() == 0)
        {
            nCells = 0;
        }
        else
        {
            labelIOList neighbour(IOobject(meshIO, "neighbour"));
            nCells = max(max(owner), max(neighbour))+1;
        }

        polyBoundaryMeshEntries patchEntries(IOobject(meshIO, "boundary"));

        basePatchNames.setSize(patchEntries.size());
        basePatchTypes.setSize(patchEntries.size());
        basePatchSizes.setSize(patchEntries.size());
        basePatchStarts.setSize(patchEntries.size());
        basePatchGroups.setSize(patchEntries.size());
        forAll(patchEntries, patchi)
        {
            basePatchNames[patchi] = patchEntries[patchi].keyword();
            const dictionary& d = patchEntries[patchi].dict();
            d.lookup("type") >> basePatchTypes[patchi];
            basePatchSizes[patchi] = readLabel(d.lookup("nFaces"));
            basePatchStarts[patchi] = readLabel(d.lookup("startFace"));
            d.readIfPresent("inGroups", basePatchGroups[patchi]);
        }
    }

    // Construct the (unallocated) mesh
    return autoPtr<unallocatedFvMesh>
    (
        new unallocatedFvMesh
        (
            io.name(),
            io.db(),                    // database
            nCells,                     // nCells
            basePatchNames,
            basePatchTypes,
            basePatchSizes,
            basePatchStarts,
            basePatchGroups,
            *reinterpret_cast<const globalMeshData*>(0)
        )
    );
}


// ************************************************************************* //
