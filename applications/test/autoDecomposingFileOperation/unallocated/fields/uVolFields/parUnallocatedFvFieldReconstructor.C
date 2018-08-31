/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "parUnallocatedFvFieldReconstructor.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::parUnallocatedFvFieldReconstructor::createPatchFaceMaps()
{
    const unallocatedFvBoundaryMesh& fvb = procMesh_.boundary();

//     // 1. Swap global cells for processor patches
//     const globalIndex globalCells(procMesh_.nCells());
//
//     labelListList remoteGlobalCells;
//     {
//         labelListList myGlobalCells(Pstream::nProcs());
//         forAll(fvb, patchI)
//         {
//             const unallocatedGenericFvPatch& pp = fvb[patchI];
//             if (pp.type() == processorFvPatch::typeName)
//             {
//                 // Use the dictionary to lookup info. Saves having full
//                 // virtual mechanism ...
//                 label nbrProci = readLabel(pp.dict().lookup("neighbProcNo"));
//
//                 const labelUList& fc = pp.faceCells();
//                 labelList& globalFc = myGlobalCells[nbrProci];
//                 globalFc.setSize(fc.size());
//                 forAll(fc, i)
//                 {
//                     globalFc[i] = globalCells.toGlobal(fc[i]);
//                 }
//             }
//         }
//
//         labelList myGlobalSizes(myGlobalCells.size());
//         forAll(myGlobalCells, proci)
//         {
//             myGlobalSizes[proci] = myGlobalCells[proci].size();
//         }
//
//         Pstream::exchange<labelList, label>
//         (
//             myGlobalCells,
//             myGlobalSizes,
//             remoteGlobalCells
//         );
//     }
//
//     // Push master faces to proc meshes
//     labelList usedMasterFaces(identity(baseMesh_.nFaces()));
//
//     // Get the master faces I am using
//distMap_.faceMap().reverseDistribute(procMesh_.nFaces(), usedMasterFaces);
//
//     // Find out to which processor each face needs to be sent. Note: not
//     // very efficient. Problem comes from single internal face becoming
//     // two processor faces on different processors. To be looked at.
//
//     const globalIndex globalProcFaces(procMesh_.nFaces());
//
//     labelList remoteMinProcFace(procMesh_.nFaces(), labelMax);
//     labelList remoteMaxProcFace(procMesh_.nFaces(), -1);
//     {
//         forAll(fvb, patchI)
//         {
//             if (fvb[patchI].type() == processorFvPatch::typeName)
//             {
//                 forAll(fvb[patchI], i)
//                 {
//                     label facei = fvb[patchI].start()+i;
//                     label globalFacei = globalProcFaces.toGlobal(facei);
//                     remoteMinProcFace[facei] = globalFacei;
//                     remoteMaxProcFace[facei] = globalFacei;
//                 }
//             }
//         }
//         mapDistributeBase::distribute
//         (
//             Pstream::commsTypes::nonBlocking,
//             List<labelPair>(),
//             distMap_.faceMap().constructSize(),
//             distMap_.faceMap().subMap(),
//             distMap_.faceMap().subHasFlip(),
//             distMap_.faceMap().constructMap(),
//             distMap_.faceMap().constructHasFlip(),
//             remoteMinProcFace,
//             minEqOp<label>(),
//             flipOp(),
//             labelMax
//         );
//         mapDistributeBase::distribute
//         (
//             Pstream::commsTypes::nonBlocking,
//             List<labelPair>(),
//             distMap_.faceMap().constructSize(),
//             distMap_.faceMap().subMap(),
//             distMap_.faceMap().subHasFlip(),
//             distMap_.faceMap().constructMap(),
//             distMap_.faceMap().constructHasFlip(),
//             remoteMaxProcFace,
//             maxEqOp<label>(),
//             flipOp(),
//             -1
//         );
//     }
//
//     labelListList sendFaces(Pstream::nProcs());
//     for (label proci = 0; proci < Pstream::nProcs(); proci++)
//     {
//         DynamicList<label> sendToProc;
//         forAll(remoteMinProcFace, baseFacei)
//         {
//             label procFacei = remoteMinProcFace[baseFacei];
//             if (globalProcFaces.isLocal(proci, procFacei))
//             {
//                 Pout<< "Found min face:"
//                     << globalProcFaces.toLocal(proci, procFacei)
//                     << " from processor:" << proci << endl;
//                 sendToProc.append(baseFacei);
//             }
//         }
//         forAll(remoteMaxProcFace, baseFacei)
//         {
//             label procFacei = remoteMaxProcFace[baseFacei];
//             if (globalProcFaces.isLocal(proci, procFacei))
//             {
//                 Pout<< "Found max face:"
//                     << globalProcFaces.toLocal(proci, procFacei)
//                     << " from processor:" << proci << endl;
//                 sendToProc.append(baseFacei);
//             }
//         }
//         sendFaces[proci].transfer(sendToProc);
//     }



//     labelListList remoteProcPatch(Pstream::nProcs());
//     forAll(remoteProcPatch, proci)
//     {
//         labelList& procPatchFaces = remoteProcPatch[proci];
//
//         // Default: return -1
//         procPatchFaces.setSize(procMesh_.nFaces(), -1);
//
//         // But on my processor return global face number
//         if (proci == Pstream::myProcNo())
//         {
//             forAll(fvb, patchI)
//             {
//                 if (fvb[patchI].type() == processorFvPatch::typeName)
//                 {
//                     forAll(fvb, i)
//                     {
//                         label facei = fvb.start()+i;
//                         procPatchFaces[facei] =
//                             globalProcFaces.toGlobal(facei);
//                     }
//                 }
//             }
//         }
//         distMap_.faceMap().distribute(procPatchFaces);
//     }
//


    // 2. Build face maps :
    //      - normal patches: boundary face to boundary face
    //      - proc patches  : for reconstruction: decomposed back to original
    patchReconFaceMaps_.setSize(baseMesh_.boundary().size());
    patchDecompFaceMaps_.setSize(baseMesh_.boundary().size());
    forAll(fvb, patchI)
    {
        if (patchI < baseMesh_.boundary().size())
        {
            // Create map for patch faces only

            // Mark all used elements (i.e. destination patch faces)
            boolList faceIsUsed(distMap_.faceMap().constructSize(), false);

            const unallocatedGenericFvPatch& basePatch =
                baseMesh_.boundary()[patchI];

            Pout<< "reconstruct map patch:" << basePatch.name()
                << " constructSize:" << distMap_.faceMap().constructSize()
                << " patch size:" << basePatch.size()
                << " patch start:" << basePatch.start()
                << endl;

            forAll(basePatch, i)
            {
                faceIsUsed[basePatch.start()+i] = true;
            }

            // Copy face map
            patchReconFaceMaps_.set
            (
                patchI,
                new mapDistributeBase(distMap_.faceMap())
            );

            // Compact out unused elements
            mapDistributeBase& map = patchReconFaceMaps_[patchI];
            labelList oldToNewSub;
            labelList oldToNewConstruct;
            map.compact
            (
                faceIsUsed,
                procMesh_.nFaces(),      // maximum index of subMap
                oldToNewSub,
                oldToNewConstruct,
                UPstream::msgType()
            );

            // Create reverse map : from undecomposed patch to proc patch
            patchDecompFaceMaps_.set
            (
                patchI,
                new mapDistributeBase
                (
                    procMesh_.boundary()[patchI].size(),
                    xferCopy(map.constructMap()),
                    xferCopy(map.subMap()),
                    map.constructHasFlip(),
                    map.subHasFlip()
                )
            );
        }
//         else if (fvb[patchI].type() == processorFvPatch::typeName)
//         {
//             // Mark all used elements (i.e. destination patch faces)
//             boolList faceIsUsed(distMap_.faceMap().constructSize(), false);
//
//             const unallocatedGenericFvPatch& procPatch =
//                 procMesh_.boundary()[patchI];
//
//             Pout<< "reconstruct map patch:" << procPatch.name()
//                 << " constructSize:" << distMap_.faceMap().constructSize()
//                 << " patch size:" << procPatch.size()
//                 << " patch start:" << procPatch.start()
//                 << endl;
//             forAll(procPatch, i)
//             {
//                 faceIsUsed[procPatch.start()+i] = true;
//             }
//
//             // Copy face map
//             patchDecompFaceMaps_.set
//             (
//                 patchI,
//                 new mapDistributeBase(distMap_.faceMap())
//             );
//
//             // Compact out unused elements
//             labelList oldToNewSub;
//             labelList oldToNewConstruct;
//             patchDecompFaceMaps_[patchI].compact
//             (
//                 faceIsUsed,
//                 procMesh_.nFaces(),      // maximum index of subMap
//                 oldToNewSub,
//                 oldToNewConstruct,
//                 UPstream::msgType()
//             );
//             //Pout<< "PROC patch:" << procPatch.name()
//             //    << " patchMap:" << patchDecompFaceMaps_[patchI] << endl;
//         }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parUnallocatedFvFieldReconstructor::parUnallocatedFvFieldReconstructor
(
    const unallocatedFvMesh& baseMesh,
    const unallocatedFvMesh& procMesh,
    const mapDistributePolyMesh& distMap
)
:
    baseMesh_(baseMesh),
    procMesh_(procMesh),
    distMap_(distMap)
{
    createPatchFaceMaps();
}


// ************************************************************************* //
