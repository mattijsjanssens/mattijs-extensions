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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::parUnallocatedFvFieldReconstructor::createPatchFaceMaps()
{
    const fvBoundaryMesh& fvb = procMesh_.boundary();

    patchFaceMaps_.setSize(fvb.size());
    forAll(fvb, patchI)
    {
        if (!isA<processorFvPatch>(fvb[patchI]))
        {
            // Create map for patch faces only

            // Mark all used elements (i.e. destination patch faces)
            boolList faceIsUsed(distMap_.faceMap().constructSize(), false);

            const unallocatedGenericFvPatch& basePatch =
                baseMesh_.boundary()[patchI];

            forAll(basePatch, i)
            {
                faceIsUsed[basePatch.start()+i] = true;
            }

            // Copy face map
            patchFaceMaps_.set
            (
                patchI,
                new mapDistributeBase(distMap_.faceMap())
            );

            // Compact out unused elements
            labelList oldToNewSub;
            labelList oldToNewConstruct;
            patchFaceMaps_[patchI].compact
            (
                faceIsUsed,
                procMesh_.nFaces(),      // maximum index of subMap
                oldToNewSub,
                oldToNewConstruct,
                UPstream::msgType()
            );
            //Pout<< "patchMap:" << patchFaceMaps_[patchI] << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parUnallocatedFvFieldReconstructor::parUnallocatedFvFieldReconstructor
(
    unallocatedFvMesh& baseMesh,
    const fvMesh& procMesh,
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
