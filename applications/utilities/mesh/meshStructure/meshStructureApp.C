/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 M Janssens
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
    meshStructure

Description
    Outputs the ordering of the mesh (if any)

\*---------------------------------------------------------------------------*/

#include "argList.H"
//#include "IOobjectList.H"
#include "fvMesh.H"
// #include "polyTopoChange.H"
// #include "ReadFields.H"
#include "volFields.H"
// #include "surfaceFields.H"
// #include "SortableList.H"
// #include "decompositionMethod.H"
// #include "renumberMethod.H"
#include "zeroGradientFvPatchFields.H"
// #include "CuthillMcKeeRenumber.H"
#include "meshStructure.H"
// #include "cellSet.H"
// #include "faceSet.H"
// #include "pointSet.H"
#include "uindirectPrimitivePatch.H"
#include "topoDistanceData.H"
#include "OppositeFaceCellWave.H"

using namespace Foam;


// Create named field from labelList for postprocessing
tmp<volScalarField> createScalarField
(
    const fvMesh& mesh,
    const word& name,
    const labelList& elems
)
{
    tmp<volScalarField> tfld
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& fld = tfld.ref();

    forAll(fld, celli)
    {
       fld[celli] = elems[celli];
    }

    return tfld;
}


void opposingFaceLabels
(
    const primitiveMesh& mesh,
    const label cellI,
    const label masterFaceLabel,
    DynamicList<label>& oppositeFaceLabels
)
{
    // Variant of cell::opposingFaceLabel

    // Algorithm:
    // Go through all the faces of the cell and find the one which
    // does not share a single vertex with the master face.  If there
    // are two or more such faces, return the first one and issue a
    // warning; if there is no opposite face, return -1;

    const face& masterFace = mesh.faces()[masterFaceLabel];

    const labelList& curFaceLabels = mesh.cells()[cellI];

    oppositeFaceLabels.clear();

    forAll(curFaceLabels, facei)
    {
        // Compare the face with the master
        const face& curFace = mesh.faces()[curFaceLabels[facei]];

        // Skip the master face
        if (curFaceLabels[facei] != masterFaceLabel)
        {
            bool sharedPoint = false;

            // Compare every vertex of the current face agains the
            // vertices of the master face
            forAll(curFace, pointi)
            {
                const label l = curFace[pointi];

                forAll(masterFace, masterPointi)
                {
                    if (masterFace[masterPointi] == l)
                    {
                        sharedPoint = true;
                        break;
                    }
                }

                if (sharedPoint) break;
            }

            // If no points are shared, this is the opposite face
            if (!sharedPoint)
            {
                // Found opposite face
                oppositeFaceLabels.append(curFaceLabels[facei]);
            }
        }
    }
}


// void handleProcPatches
// (
//     const polyMesh& mesh,
//     const List<topoDistanceData>& allFaceInfo,
//     DynamicList<label>& frontFaces
// )
// {
//     const globalMeshData& pData = mesh.globalData();
// 
//     // Which patches are processor patches
//     const labelList& procPatches = pData.processorPatches();
// 
//     PackedBoolList isFrontFace(mesh.nFaces()-mesh.nInternalFaces());
//     forAll(frontFaces, i)
//     {
//         if (!mesh.isInternalFace(frontFaces[i]))
//         {
//             isFrontFace[frontFaces[i]-mesh.nInternalFaces()] = true;
//         }
//     }
// 
// 
//     // Send all
// 
//     PstreamBuffers pBufs(Pstream::nonBlocking);
// 
//     forAll(procPatches, i)
//     {
//         label patchi = procPatches[i];
// 
//         const processorPolyPatch& procPatch =
//             refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);
// 
//         // Allocate buffers
//         DynamicList<label> sendFaces(procPatch.size());
//         DynamicList<Type> sendFacesInfo(procPatch.size());
// 
//         // Determine which faces changed on current patch
//         forAll(procPatch, patchFacei)
//         {
//             label meshFacei = procPatch.start()+patchFacei;
// 
//             if (isFrontFace[meshFacei-mesh.nInternalFaces()])
//             {
//                 sendFaces.append(patchFacei);   // patch face
//                 sendFacesInfo.append(allFaceInfo[meshFacei]);
//             }
//         }
// 
// 
//         nSendFaces = getChangedPatchFaces
//         (
//             procPatch,
//             0,
//             procPatch.size(),
//             sendFaces,
//             sendFacesInfo
//         );
// 
//         // Adapt wallInfo for leaving domain
//         leaveDomain
//         (
//             procPatch,
//             nSendFaces,
//             sendFaces,
//             sendFacesInfo
//         );
// 
//         if (debug & 2)
//         {
//             Pout<< " Processor patch " << patchi << ' ' << procPatch.name()
//                 << " communicating with " << procPatch.neighbProcNo()
//                 << "  Sending:" << nSendFaces
//                 << endl;
//         }
// 
//         UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
//         //writeFaces(nSendFaces, sendFaces, sendFacesInfo, toNeighbour);
//         toNeighbour
//             << SubList<label>(sendFaces, nSendFaces)
//             << SubList<Type>(sendFacesInfo, nSendFaces);
//     }
// 
//     pBufs.finishedSends();
// 
//     // Receive all
// 
//     forAll(procPatches, i)
//     {
//         label patchi = procPatches[i];
// 
//         const processorPolyPatch& procPatch =
//             refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchi]);
// 
//         // Allocate buffers
//         labelList receiveFaces;
//         List<Type> receiveFacesInfo;
// 
//         {
//             UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);
//             fromNeighbour >> receiveFaces >> receiveFacesInfo;
//         }
// 
//         if (debug & 2)
//         {
//             Pout<< " Processor patch " << patchi << ' ' << procPatch.name()
//                 << " communicating with " << procPatch.neighbProcNo()
//                 << "  Receiving:" << receiveFaces.size()
//                 << endl;
//         }
// 
//         // Apply transform to received data for non-parallel planes
//         if (!procPatch.parallel())
//         {
//             transform
//             (
//                 procPatch.forwardT(),
//                 receiveFaces.size(),
//                 receiveFacesInfo
//             );
//         }
// 
//         // Adapt wallInfo for entering domain
//         enterDomain
//         (
//             procPatch,
//             receiveFaces.size(),
//             receiveFaces,
//             receiveFacesInfo
//         );
// 
//         // Merge received info
//         mergeFaceInfo
//         (
//             procPatch,
//             receiveFaces.size(),
//             receiveFaces,
//             receiveFacesInfo
//         );
//     }
// }
// 
// 
// void syncFront
// (
//     const polyMesh& mesh,
//     const List<topoDistanceData>& allFaceInfo,
//     DynamicList<label>& frontFaces
// )
// {
//     // Modify front to include coupled faces
//     
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Determines mesh structure");
    argList::validArgs.append("patches");

    #include "addRegionOption.H"
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    // Get times list
    instantList Times = runTime.times();

    // Set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    #include "createNamedMesh.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Find set of patches from the list of regular expressions provided
    const wordReList patches((IStringStream(args[1])()));
    const labelHashSet patchSet(pbm.patchSet(patches));

    label nFaces = 0;
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        nFaces += pbm[iter.key()].size();
    }

    labelList meshFaces(nFaces);
    nFaces = 0;
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        const polyPatch& pp = pbm[iter.key()];
        forAll(pp, i)
        {
            meshFaces[nFaces++] = pp.start()+i;
        }
    }

    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(mesh.faces(), meshFaces),
        mesh.points()
    );

Pout<< "pp:" << pp.size() << endl;

    //meshStructure ms(mesh, upp);

    List<topoDistanceData> allCellInfo(mesh.nCells());
    List<topoDistanceData> allFaceInfo(mesh.nFaces());
//     {
//         // Start of changes
//         DynamicList<label> frontFaces(pp.size());
//         DynamicList<topoDistanceData> frontData(pp.size());
//         forAll(pp, patchFacei)
//         {
//             label meshFacei = pp.addressing()[patchFacei];
// 
//             // Make sure face present only once in initial front
//             if (!allFaceInfo[meshFacei].valid(dummyTrackData))
//             {
//                 allFaceInfo[meshFacei] = topoDistanceData(patchFacei, 0);
//                 frontData.append(allFaceInfo[meshFacei]);
//                 frontFaces.append(meshFacei);
//             }
//         }
// 
// 
//         PackedBoolList isNewFrontFace(mesh.nFaces());
// 
//         DynamicList<label> oppositeFaceLabels;
// 
//         while (true)
//         {
//             DynamicList<label> newFrontFaces(frontFaces.size());
//             isNewFrontFace = false;
// 
//             // Collect opposite face
//             forAll(frontFaces, i)
//             {
//                 label facei = frontFaces[i];
//                 Pout<< "face:" << facei << " at:" << mesh.faceCentres()[facei]
//                     << " distance:" << allFaceInfo[facei].distance()
//                     << endl;
// 
//                 {
//                     // Owner side
// 
//                     label own = mesh.faceOwner()[facei];
//                     opposingFaceLabels
//                     (
//                         mesh,
//                         own,
//                         facei,
//                         oppositeFaceLabels
//                     );
//                     Pout<< "    own:" << own
//                         << " at:" << mesh.cellCentres()[own]
//                         << " oppositefaces:"
//                         << pointField(mesh.faceCentres(), oppositeFaceLabels)
//                         << endl;
// 
//                     if (oppositeFaceLabels.size())
//                     {
//                         bool propagate = allCellInfo[own].updateCell
//                         (
//                             mesh,
//                             own,
//                             facei,
//                             allFaceInfo[facei],
//                             1e-6,                   //tol,
//                             dummyTrackData
//                         );
// 
//                         if (propagate && oppositeFaceLabels.size() == 1)
//                         {
//                             label oppFacei = oppositeFaceLabels[0];
// 
//                             bool propagate = allFaceInfo[oppFacei].updateFace
//                             (
//                                 mesh,
//                                 oppFacei,
//                                 own,
//                                 allCellInfo[own],
//                                 1e-6,           //tol,
//                                 dummyTrackData
//                             );
//                             if (propagate && isNewFrontFace.set(oppFacei))
//                             {
//                                 newFrontFaces.append(oppFacei);
//                             }
//                         }
//                     }
//                 }
//                 if (mesh.isInternalFace(facei))
//                 {
//                     label nei = mesh.faceNeighbour()[facei];
//                     opposingFaceLabels
//                     (
//                         mesh,
//                         nei,
//                         facei,
//                         oppositeFaceLabels
//                     );
//                     Pout<< "    nei:" << nei
//                         << " at:" << mesh.cellCentres()[nei]
//                         << " oppositefaces:"
//                         << pointField(mesh.faceCentres(), oppositeFaceLabels)
//                         << endl;
// 
// 
//                     if (oppositeFaceLabels.size())
//                     {
//                         bool propagate = allCellInfo[nei].updateCell
//                         (
//                             mesh,
//                             nei,
//                             facei,
//                             allFaceInfo[facei],
//                             1e-6,               //tol,
//                             dummyTrackData
//                         );
// 
//                         if (propagate && oppositeFaceLabels.size() == 1)
//                         {
//                             label oppFacei = oppositeFaceLabels[0];
// 
//                             bool propagate = allFaceInfo[oppFacei].updateFace
//                             (
//                                 mesh,
//                                 oppFacei,
//                                 nei,
//                                 allCellInfo[nei],
//                                 1e-6,                   //tol,
//                                 dummyTrackData
//                             );
//                             if (propagate && isNewFrontFace.set(oppFacei))
//                             {
//                                 newFrontFaces.append(oppFacei);
//                             }
//                         }
//                     }
//                 }
//             }
// 
// 
//             //TBD. Synchronise across coupled faces ...
// 
//             Pout<< "old front:" << frontFaces
//                 << " new front:" << newFrontFaces << endl;
// 
//             if (returnReduce(newFrontFaces.size(), sumOp<label>()) == 0)
//             {
//                 break;
//             }
// 
//             frontFaces.transfer(newFrontFaces);
// 
//             // Do parallel synchronisation
//             syncFront(mesh, allFaceInfo, frontFaces);
//         }
//     }



    DynamicList<label> patchFaces(pp.size());
    DynamicList<topoDistanceData> patchData(pp.size());
    forAll(pp, patchFacei)
    {
        patchFaces.append(pp.addressing()[patchFacei]);
        patchData.append(topoDistanceData(patchFacei, 0));
    }

Pout<< " patchFaces:" << patchFaces << endl;

    OppositeFaceCellWave<topoDistanceData> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        allFaceInfo,
        allCellInfo,
        mesh.globalData().nTotalCells()+1
    );

    volScalarField fld
    (
        IOobject
        (
            "distance",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    forAll(fld, celli)
    {
       fld[celli] = allCellInfo[celli].distance();
    }
    fld.correctBoundaryConditions();
    fld.write();


//     createScalarField
//     (
//         mesh,
//         "cellLayer",
//         ms.cellLayer()
//     )().write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
