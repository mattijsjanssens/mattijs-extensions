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

#include "dynamicPatchSplitFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "sigFpe.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicPatchSplitFvMesh, 0);
    addToRunTimeSelectionTable
    (dynamicFvMesh, dynamicPatchSplitFvMesh, IOobject);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicPatchSplitFvMesh::dynamicPatchSplitFvMesh(const IOobject& io)
:
    dynamicFvMesh(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicPatchSplitFvMesh::~dynamicPatchSplitFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicPatchSplitFvMesh::update()
{
    const label patchi = 2;
    static label patchFacei = 0;

    const polyPatch& pp = boundaryMesh()[patchi];

    boolList isBoundarySplit(nFaces()-nInternalFaces(), false);

    // Insert face to split
    isBoundarySplit[pp.start()+patchFacei-nInternalFaces()] = true;

    // Make sure any coupled faces are also split
    syncTools::syncBoundaryFaceList(*this, isBoundarySplit, orEqOp<bool>());


    polyTopoChange meshMod(*this);

    forAll(isBoundarySplit, bFacei)
    {
        if (isBoundarySplit[bFacei])
        {
            label meshFacei = nInternalFaces()+bFacei;

            Pout<< "Splitting face:" << meshFacei
                << " patch:" << boundaryMesh().whichPatch(meshFacei)
                << " at " << faceCentres()[meshFacei] << endl;

            const face& f = faces()[meshFacei];

            // Add face centre
            label pointi = meshMod.addPoint
            (
                faceCentres()[meshFacei],
                -1,
                -1,
                true
            );

            // Remove existing face
            meshMod.removeFace(meshFacei, -1);

            // Add triangles
            forAll(f, fp)
            {
                face tri(3);
                tri[0] = f[fp];
                tri[1] = f[f.fcIndex(fp)];
                tri[2] = pointi;
                meshMod.addFace
                (
                    tri,
                    faceOwner()[meshFacei],
                    -1,                         // nei
                    -1,                         // masterPointID
                    -1,                         // masterEdgeID
                    meshFacei,                  // masterFaceID
                    false,                      // flip
                    boundaryMesh().whichPatch(meshFacei),
                    -1,                         // zone
                    false                       // zone sign
                );

                Pout<< "Added face:" << tri
                    << " verts:"
                    << UIndirectList<point>(meshMod.points(), tri)
                    << endl;
            }
        }
    }

    // Change mesh and generate map.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);

    // Update fields
    updateMesh(map);

    // Move mesh
    if (map().hasMotionPoints())
    {
// Note: added points:
//  - mesh.points() = vector::zero
//  - preMotionPoints = wanted location (face centre)
// Tbd. have an 'old' position in the call to meshMod.addPoint.
DebugVar(map().preMotionPoints()-points());

        movePoints(map().preMotionPoints());
    }

    patchFacei++;


    const bool hasChanged = true;

    topoChanging(hasChanged);
    moving(hasChanged);

    return hasChanged;
}


// ************************************************************************* //
