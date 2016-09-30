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
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(smoothTriSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, smoothTriSurfaceMesh, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void Foam::smoothTriSurfaceMesh::calcPointNormals
// (
//     const PackedBoolList& isBorderEdge,
//     vectorField& pointNormals
// )
// {
//     const triSurface& s = static_cast<const triSurface&>(*this);
//     const pointField& pts = s.points();
// 
//     // Calculate faceNormals. Could use cached version from underlying
//     // PrimitivePatch but want to avoid extra storage.
//     vectorField faceNormals(size());
//     forAll(faceNormals, facei)
//     {
//         faceNormals[facei] = s[facei].normal(pts)
//         faceNormals[i] /= mag(faceNormals[i]) + VSMALL;
//     }
// 
//     const labelListList& faceEdges = s.faceEdges();
//     const labelListList& edgeFaces = s.edgeFaces();
//     const edgeList& edges = s.edges();
// 
// 
// 
//     // Calculate average edge normals
//     vectorField edgeNormals(pp.nEdges(), Zero);
//     forAll(edgeNormals, edgei)
//     {
//         const labelList& eFaces = edgeFaces[edgei];
//         forAll(eFaces, i)
//         {
//             edgeNormals[edgei] += faceNormals[eFaces[i]];
//         }
//     }
//     edgeNormals /= mag(edgeNormals)+VSMALL;
// 
// 
// 
//     pointNormals.setSize(s.nPoints());
//     pointNormals = Zero;
// 
//     forAll(faceEdges, facei)
//     {
//         const labelList& fEdges = faceEdges[facei];
// 
//         forAll(fEdges, i)
//         {
//             label edgei = fEdges[i];
//             vector edgeN;
//             if (isBorderEdge[edgei])
//             {
//                 edgeN = faceNormals[facei];
//             }
//             else
//             {
//                 edgeN = edgeNormals[edgei];
//             }
// 
//             const edge& e = edges[edgei];
//             pointNormals[e[0]] += edgeN;
//             pointNormals[e[1]] += edgeN;
//         }
//     }
// 
//     pointNormals /= mag(pointNormals)+VSMALL;
// }

void Foam::smoothTriSurfaceMesh::calcFeatureEdges
(
    const scalar featureAngle
)
{
    scalar cosAngle = Foam::cos(degToRad(featureAngle));

    const triSurface& s = static_cast<const triSurface&>(*this);
    const labelListList& edgeFaces = s.edgeFaces();
    const vectorField& faceNormals = s.faceNormals();

    forAll(edgeFaces, edgei)
    {
        const labelList& eFaces = edgeFaces[edgei];

        if (eFaces.size() > 2)
        {
            isBorderEdge_[edgei] = true;
        }
        else if (eFaces.size() == 2)
        {
            const vector& n0 = faceNormals[eFaces[0]];
            const vector& n1 = faceNormals[eFaces[1]];
            if ((n0&n1) < cosAngle)
            {
                isBorderEdge_[edgei] = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
// (
//     const IOobject& io,
//     const triSurface& s
// )
// :
//     triSurfaceMesh(io, s)
// {}


Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
(
    const IOobject& io,
    const scalar featureAngle
)
:
    triSurfaceMesh(io),
    isBorderEdge_(nEdges())
{
    calcFeatureEdges(featureAngle);
}


Foam::smoothTriSurfaceMesh::smoothTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, dict),
    isBorderEdge_(nEdges())
{
    if (dict.found("featureAngle"))
    {
        calcFeatureEdges(readScalar(dict.lookup("featureAngle")));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smoothTriSurfaceMesh::~smoothTriSurfaceMesh()
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
    const labelListList& faceEdges = s.faceEdges();
    const labelListList& edgeFaces = s.edgeFaces();
    const vectorField& faceNormals = s.faceNormals();

    normal.setSize(info.size());

    forAll(info, i)
    {
        if (info[i].hit())
        {
            label facei = info[i].index();
            const labelList& fEdges = faceEdges[facei];

            // Get local point normals
            FixedList<vector, 3> n(Zero);
            forAll(fEdges, fEdgei)
            {
                label edgei = fEdges[fEdgei];

                vector edgeN;
                if (!isBorderEdge_[edgei])
                {
                    const labelList& eFaces = edgeFaces[edgei];

                    edgeN = Zero;
                    forAll(eFaces, eFacei)
                    {
                        edgeN += faceNormals[eFaces[eFacei]];
                    }
                    edgeN /= mag(edgeN) + VSMALL;
                }
                else
                {
                    edgeN = faceNormals[facei];
                }
                n[fEdgei] += edgeN;
                n[n.fcIndex(fEdgei)] += edgeN;
            }

            forAll(n, ni)
            {
                n[ni] /= mag(n[ni]) + VSMALL;
            }

            // Get local coordinates in triangle
            FixedList<scalar, 3> coordinates;
            s[facei].tri(pts).barycentric(info[i].hitPoint(), coordinates);

            forAll(coordinates, ci)
            {
                normal[i] = coordinates[ci]*n[ci];
                normal[i] /= mag(normal[i]) + VSMALL;
            }

            Pout<< "face:" << facei << nl
                << "    fc:" << faceCentres()[facei] << nl
                << "    n:" << faceNormals[facei] << nl
                << "    smooth:" << normal[i] << nl
                << endl;
        }
        else
        {
            // Set to what?
            normal[i] = Zero;
        }
    }
}


// ************************************************************************* //
