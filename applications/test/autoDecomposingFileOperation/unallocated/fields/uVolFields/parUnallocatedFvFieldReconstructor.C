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
    // Build face maps : decomposed boundary face to undecomposed boundary
    //                   face.
    // Note that processor patches get handled differently
    patchReconFaceMaps_.setSize(baseMesh_.boundary().size());
    patchDecompFaceMaps_.setSize(baseMesh_.boundary().size());
    forAll(baseMesh_.boundary(), patchI)
    {
        // Mark all used elements (i.e. destination patch faces)
        boolList faceIsUsed(distMap_.faceMap().constructSize(), false);

        const unallocatedGenericFvPatch& basePatch =
            baseMesh_.boundary()[patchI];

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
