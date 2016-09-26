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

\*---------------------------------------------------------------------------*/

#include "smoothTriSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "EdgeMap.H"
#include "triSurfaceFields.H"
#include "Time.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(smoothTriSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, smoothTriSurfaceMesh, dict);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
(
    const IOobject& io,
    const triSurface& s
)
:
    triSurfaceMesh(io, s)
{}


Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh(const IOobject& io)
:
    triSurfaceMesh(io)
{}


Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smoothTriSurfaceMesh::~triSurfaceMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoothTriSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    const triSurface& s = static_cast<const triSurface&>(*this);
    const pointField& pts = s.points();

    normal.setSize(info.size());

    if (minQuality_ >= 0)
    {
        // Make sure we don't use triangles with low quality since
        // normal is not reliable.

        const labelListList& faceFaces = s.faceFaces();

        forAll(info, i)
        {
            if (info[i].hit())
            {
                label facei = info[i].index();
                normal[i] = s[facei].normal(pts);

                scalar qual = s[facei].tri(pts).quality();

                if (qual < minQuality_)
                {
                    // Search neighbouring triangles
                    const labelList& fFaces = faceFaces[facei];

                    forAll(fFaces, j)
                    {
                        label nbrI = fFaces[j];
                        scalar nbrQual = s[nbrI].tri(pts).quality();
                        if (nbrQual > qual)
                        {
                            qual = nbrQual;
                            normal[i] = s[nbrI].normal(pts);
                        }
                    }
                }

                normal[i] /= mag(normal[i]) + VSMALL;
            }
            else
            {
                // Set to what?
                normal[i] = Zero;
            }
        }
    }
    else
    {
        forAll(info, i)
        {
            if (info[i].hit())
            {
                label facei = info[i].index();
                // Cached:
                //normal[i] = faceNormals()[facei];

                // Uncached
                normal[i] = s[facei].normal(pts);
                normal[i] /= mag(normal[i]) + VSMALL;
            }
            else
            {
                // Set to what?
                normal[i] = Zero;
            }
        }
    }
}


void Foam::smoothTriSurfaceMesh::setField(const labelList& values)
{
    autoPtr<triSurfaceLabelField> fldPtr
    (
        new triSurfaceLabelField
        (
            IOobject
            (
                "values",
                objectRegistry::time().timeName(),  // instance
                "triSurface",                       // local
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimless,
            labelField(values)
        )
    );

    // Store field on triMesh
    fldPtr.ptr()->store();
}


void Foam::smoothTriSurfaceMesh::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    if (foundObject<triSurfaceLabelField>("values"))
    {
        values.setSize(info.size());

        const triSurfaceLabelField& fld = lookupObject<triSurfaceLabelField>
        (
            "values"
        );

        forAll(info, i)
        {
            if (info[i].hit())
            {
                values[i] = fld[info[i].index()];
            }
        }
    }
}


void Foam::smoothTriSurfaceMesh::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        if (!tree().bb().contains(pt))
        {
            // Have to calculate directly as outside the octree
            volType[pointi] = tree().shapes().getVolumeType(tree(), pt);
        }
        else
        {
            // - use cached volume type per each tree node
            volType[pointi] = tree().getVolumeType(pt);
        }
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


bool Foam::smoothTriSurfaceMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    fileName fullPath(searchableSurface::objectPath());

    if (!mkDir(fullPath.path()))
    {
        return false;
    }

    triSurface::write(fullPath);

    if (!isFile(fullPath))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //
