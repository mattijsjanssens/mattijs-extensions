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

#include "isoSurfaceTopo.H"
#include "polyMesh.H"
#include "tetMatcher.H"
#include "tetPointRef.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoSurfaceTopo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::isoSurfaceTopo::isTriCut
(
    const triFace& tri,
    const scalarField& pointValues
) const
{
    bool aLower = (pointValues[tri[0]] < iso_);
    bool bLower = (pointValues[tri[1]] < iso_);
    bool cLower = (pointValues[tri[2]] < iso_);

    return !(aLower == bLower && aLower == cLower);
}


Foam::isoSurfaceTopo::cellCutType Foam::isoSurfaceTopo::calcCutType
(
    const bool isTet,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];

    if (isTet)
    {
        forAll(cFaces, cFacei)
        {
            const face& f = mesh_.faces()[cFaces[cFacei]];

            for (label fp = 1; fp < f.size() - 1; fp++)
            {
                triFace tri(f[0], f[fp], f[f.fcIndex(fp)]);

                if (isTriCut(tri, pVals_))
                {
                    return CUT;
                }
            }
        }
        return NOTCUT;
    }
    else
    {
        bool cellLower = (cVals_[celli] < iso_);

        // First check if there is any cut in cell
        bool edgeCut = false;

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const face& f = mesh_.faces()[facei];

            // Check pyramids cut
            forAll(f, fp)
            {
                if ((pVals_[f[fp]] < iso_) != cellLower)
                {
                    edgeCut = true;
                    break;
                }
            }

            if (edgeCut)
            {
                break;
            }

            const label fp0 = mesh_.tetBasePtIs()[facei];
            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                label nextFp = f.fcIndex(fp);

                if (isTriCut(triFace(f[fp0], f[fp], f[nextFp]), pVals_))
                {
                    edgeCut = true;
                    break;
                }

                fp = nextFp;
            }

            if (edgeCut)
            {
                break;
            }
        }

        if (edgeCut)
        {
            // Count actual cuts (expensive since addressing needed)
            // Note: not needed if you don't want to preserve maxima/minima
            // centred around cellcentre. In that case just always return CUT

            const labelList& cPoints = mesh_.cellPoints(celli);

            label nPyrCuts = 0;

            forAll(cPoints, i)
            {
                if ((pVals_[cPoints[i]] < iso_) != cellLower)
                {
                    nPyrCuts++;
                }
            }

            if (nPyrCuts == cPoints.size())
            {
                return SPHERE;
            }
            else
            {
                return CUT;
            }
        }
        else
        {
            return NOTCUT;
        }
    }
}


Foam::label Foam::isoSurfaceTopo::calcCutTypes
(
    tetMatcher& tet,
    List<cellCutType>& cellCutTypes
)
{
    const cellList& cells = mesh_.cells();

    cellCutTypes.setSize(cells.size());
    label nCutCells = 0;
    forAll(cells, celli)
    {
        cellCutTypes[celli] = calcCutType(tet.isA(mesh_, celli), celli);

        if (cellCutTypes[celli] == CUT)
        {
            nCutCells++;
        }
    }

    if (debug)
    {
        Pout<< "isoSurfaceTopo : detected " << nCutCells
            << " candidate cut cells." << endl;
    }
    return nCutCells;
}


Foam::label Foam::isoSurfaceTopo::generatePoint
(
    const label facei,
    const bool edgeIsDiag,
    const edge& vertices,

    DynamicList<edge>& pointToVerts,
    DynamicList<labelList>& pointToFace,
    DynamicList<bool>& pointFromDiag,
    EdgeMap<label>& vertsToPoint
) const
{
    EdgeMap<label>::const_iterator edgeFnd = vertsToPoint.find(vertices);
    if (edgeFnd != vertsToPoint.end())
    {
        label pointi = edgeFnd();
        labelList& pFaces = pointToFace[pointi];

//Pout<< "for vertices:" << vertices << " found already generated point "
//    << pointi << " with faces:" << pFaces << " facei:" << facei << endl;
        if (findIndex(pFaces, facei) == -1)
        {
            pFaces.append(facei);
        }
        return edgeFnd();
    }
    else
    {
        // Generate new point
        label pointi = pointToVerts.size();

        pointToVerts.append(vertices);
        pointToFace.append(labelList(1, facei));
        pointFromDiag.append(edgeIsDiag);
        vertsToPoint.insert(vertices, pointi);
        return pointi;
    }
}


void Foam::isoSurfaceTopo::generateTriPoints
(
    const label facei,
    const FixedList<scalar, 4>& s,
    const FixedList<point, 4>& p,
    const FixedList<label, 4>& pIndex,
    const FixedList<bool, 6>& edgeIsDiag,// per tet edge whether is face diag

    DynamicList<edge>& pointToVerts,
    DynamicList<labelList>& pointToFace,
    DynamicList<bool>& pointFromDiag,

    EdgeMap<label>& vertsToPoint,
    DynamicList<label>& verts       // every three verts is new triangle
) const
{
    int triIndex = 0;
    if (s[0] < iso_)
    {
        triIndex |= 1;
    }
    if (s[1] < iso_)
    {
        triIndex |= 2;
    }
    if (s[2] < iso_)
    {
        triIndex |= 4;
    }
    if (s[3] < iso_)
    {
        triIndex |= 8;
    }

    /* Form the vertices of the triangles for each case */
    switch (triIndex)
    {
        case 0x00:
        case 0x0F:
        break;

        case 0x01:
        case 0x0E:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    edge(pIndex[0], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    edge(pIndex[0], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    edge(pIndex[0], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (triIndex == 0x0E)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x02:
        case 0x0D:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    edge(pIndex[1], pIndex[0]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    edge(pIndex[1], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    edge(pIndex[1], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (triIndex == 0x0D)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x03:
        case 0x0C:
        {
            label p0p2
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    edge(pIndex[0], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            label p1p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    edge(pIndex[1], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    edge(pIndex[0], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p1p3);
            verts.append(p0p2);
            verts.append(p1p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    edge(pIndex[1], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p0p2);

            if (triIndex == 0x0C)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-5], verts[sz-4]);
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x04:
        case 0x0B:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    edge(pIndex[2], pIndex[0]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    edge(pIndex[2], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    edge(pIndex[2], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (triIndex == 0x0B)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x05:
        case 0x0A:
        {
            label p0p1
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    edge(pIndex[0], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            label p2p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    edge(pIndex[2], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append(p0p1);
            verts.append(p2p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    edge(pIndex[0], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p0p1);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    edge(pIndex[1], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p2p3);

            if (triIndex == 0x0A)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-5], verts[sz-4]);
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x06:
        case 0x09:
        {
            label p0p1
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    edge(pIndex[0], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            label p2p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    edge(pIndex[2], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append(p0p1);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    edge(pIndex[1], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p2p3);
            verts.append(p0p1);
            verts.append(p2p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    edge(pIndex[0], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (triIndex == 0x09)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-5], verts[sz-4]);
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x08:
        case 0x07:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    edge(pIndex[3], pIndex[0]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    edge(pIndex[3], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    edge(pIndex[3], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            if (triIndex == 0x07)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;
    }
}


void Foam::isoSurfaceTopo::generateTriPoints
(
    const polyMesh& mesh,
    const label celli,
    const bool isTet,

    DynamicList<edge>& pointToVerts,
    DynamicList<labelList>& pointToFace,
    DynamicList<bool>& pointFromDiag,

    EdgeMap<label>& vertsToPoint,
    DynamicList<label>& verts,
    DynamicList<label>& faceLabels
) const
{
    const cell& cFaces = mesh.cells()[celli];

    if (isTet)
    {
        // For tets don't do cell-centre decomposition, just use the
        // tet points and values

        label facei = cFaces[0];
        const face& f0 = mesh_.faces()[facei];

        // Get the other point
        const face& f1 = mesh_.faces()[cFaces[1]];
        label oppositeI = -1;
        forAll(f1, fp)
        {
            oppositeI = f1[fp];
            if (findIndex(f0, oppositeI) == -1)
            {
                break;
            }
        }


        label p0 = f0[0];
        label p1 = f0[1];
        label p2 = f0[2];
        FixedList<bool, 6> edgeIsDiag(false);

        if (mesh.faceOwner()[facei] == celli)
        {
            Swap(p1, p2);
        }

        tetPointRef tet
        (
            mesh.points()[p0],
            mesh.points()[p1],
            mesh.points()[p2],
            mesh.points()[oppositeI]
        );

        label startTrii = verts.size();
        generateTriPoints
        (
            facei,
            FixedList<scalar, 4>
            ({
                pVals_[p0],
                pVals_[p1],
                pVals_[p2],
                pVals_[oppositeI]
            }),
            FixedList<point, 4>
            ({
                mesh.points()[p0],
                mesh.points()[p1],
                mesh.points()[p2],
                mesh.points()[oppositeI]
            }),
            FixedList<label, 4>
            ({
                p0,
                p1,
                p2,
                oppositeI
            }),
            edgeIsDiag,

            pointToVerts,
            pointToFace,
            pointFromDiag,
            vertsToPoint,
            verts       // every three verts is new triangle
        );

        label nTris = (verts.size()-startTrii)/3;
        for (label i = 0; i < nTris; i++)
        {
            faceLabels.append(facei);
        }
    }
    else
    {
        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const face& f = mesh.faces()[facei];

            label fp0 = mesh.tetBasePtIs()[facei];

            label startTrii = verts.size();

            // Skip undefined tets
            if (fp0 < 0)
            {
                fp0 = 0;
            }

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                label nextFp = f.fcIndex(fp);

                FixedList<bool, 6> edgeIsDiag(false);

                label p0 = f[fp0];
                label p1 = f[fp];
                label p2 = f[nextFp];
                if (mesh.faceOwner()[facei] == celli)
                {
                    Swap(p1, p2);
                    if (i != 2) edgeIsDiag[1] = true;
                    if (i != f.size()-1) edgeIsDiag[0] = true;
                }
                else
                {
                    if (i != 2) edgeIsDiag[0] = true;
                    if (i != f.size()-1) edgeIsDiag[1] = true;
                }

                tetPointRef tet
                (
                    mesh.points()[p0],
                    mesh.points()[p1],
                    mesh.points()[p2],
                    mesh.cellCentres()[celli]
                );

                generateTriPoints
                (
                    facei,
                    FixedList<scalar, 4>
                    ({
                        pVals_[p0],
                        pVals_[p1],
                        pVals_[p2],
                        cVals_[celli]
                    }),
                    FixedList<point, 4>
                    ({
                        mesh.points()[p0],
                        mesh.points()[p1],
                        mesh.points()[p2],
                        mesh.cellCentres()[celli]
                    }),
                    FixedList<label, 4>
                    ({
                        p0,
                        p1,
                        p2,
                        mesh.nPoints()+celli
                    }),
                    edgeIsDiag,

                    pointToVerts,
                    pointToFace,
                    pointFromDiag,
                    vertsToPoint,
                    verts       // every three verts is new triangle
                );

                fp = nextFp;
            }

            label nTris = (verts.size()-startTrii)/3;
            for (label i = 0; i < nTris; i++)
            {
                faceLabels.append(facei);
            }
        }
    }
}


Foam::label Foam::isoSurfaceTopo::findMeshEdge
(
    const boolList& pointFromDiag,
    const primitivePatch& pp,
    const labelList& loop,
    const label startFp
) const
{
    const labelList& mp = pp.meshPoints();

    label fp = startFp;

    forAll(loop, i)
    {
        label pointi = mp[loop[fp]];
        if (!pointFromDiag[pointi])
        {
            return fp;
        }

        fp = loop.fcIndex(fp);
    }
    return -1;
}


void Foam::isoSurfaceTopo::detectFaceEdges
(
    const primitivePatch& pp,
    const labelList& loop,
    const boolList& pointFromDiag,
    const labelListList& pointToFace,
    EdgeMap<label>& vertsToFace
) const
{
    const labelList& mp = pp.meshPoints();

    // Search for mesh edge
    label edgeFp = findMeshEdge(pointFromDiag, pp, loop, 0);

    if (edgeFp != -1)
    {
        label fp = f.fcIndex(edgeFp);
        for (label i = 0; i < loop.size()-1; i++)
        {
            const label pointi = mp[loop[fp]];
            if (pointFromDiag[pointi])
            {
                label facei = pointToFace[pointi][0];

                // Search forwards for next meshFp
                label nextFp =
                    findMeshEdge(pointFromDiag, pp, loop, f.fcIndex(fp));

                const edge vertices
                (
                    mp[loop[edgeFp]],
                    mp[loop[nextFp]]
                );
                vertsToFaceinsert(vertices, facei);
            }
            else
            {
                edgeFp = fp;
            }

            fp = f.fcIndex(fp);
        }
    }
}


void Foam::isoSurfaceTopo::filterFace
(
    const primitivePatch& pp,
    const labelList& loop,
    const boolList& pointFromDiag,
    const labelListList& pointToFace,
    face& f,
    EdgeMap<label>& vertsToFace
) const
{
    const labelList& mp = pp.meshPoints();

    f.setSize(loop.size());
    label fpi = 0;

    // Search for mesh edge
    label edgeFp = findMeshEdge(pointFromDiag, pp, loop, 0);

    if (edgeFp == -1)
    {
        //Pout<< "**TBD." << endl;
        forAll(loop, fp)
        {
            label pointi = mp[loop[fp]];
            f[fpi++] = pointi;
        }
    }
    else
    {
        f[fpi++] = mp[loop[edgeFp]];

        label fp = f.fcIndex(edgeFp);
        for (label i = 0; i < loop.size()-1; i++)
        {
            const label pointi = mp[loop[fp]];
            if (pointFromDiag[pointi])
            {
                label facei = pointToFace[pointi][0];

                // Search forwards for next meshFp
                label nextFp =
                    findMeshEdge(pointFromDiag, pp, loop, f.fcIndex(fp));

//Pout<< "On loop:" << loop
//    << " fp:" << fp
//    << " leftFp:" << edgeFp
//    << " rightFp:" << nextFp
//    << endl;
//
//Pout<< " leftReal:" << mp[loop[edgeFp]]
//    << " at:" << pp.points()[mp[loop[edgeFp]]]
//    << " mid:" << pointi << " at:" << pp.points()[pointi]
//    << " rightReal:" << mp[loop[nextFp]]
//    << " at:" << pp.points()[mp[loop[nextFp]]]
//    << endl;

                const edge vertices
                (
                    mp[loop[edgeFp]],
                    mp[loop[nextFp]]
                );
                EdgeMap<label>::const_iterator edgeFnd =
                    vertsToFace.find(vertices);
                if (edgeFnd != vertsToFace.end())
                {
                    if (edgeFnd() != facei)
                    {
                        //Pout<< "** face:" << facei
                        //    << " would generate same edge:" << vertices
                        //    << " at:" << vertices.line(pp.points())
                        //    << " as face:" << edgeFnd() << endl;

                        // Protect i.e. keep point
                        f[fpi++] = pointi;
                    }
                }
                else
                {
                    // First occurence of vertex pair
                    vertsToFace.insert(vertices, facei);
                    //Pout<< "storing edge:" << vertices
                    //    << " at:" << vertices.line(pp.points())
                    //    << " generated from face:" << facei
                    //    << endl;
                }
            }
            else
            {
                f[fpi++] = pointi;
                edgeFp = fp;
            }

            fp = f.fcIndex(fp);
        }

        if (fpi > 2)
        {
            f.setSize(fpi);
        }
        else
        {
            // Keep original face
            forAll(f, i)
            {
                label pointi = mp[loop[i]];
                f[i] = pointi;
            }
        }
    }
}


void Foam::isoSurfaceTopo::triangulateOutside
(
    const bool filterDiag,
    const primitivePatch& pp,
    const boolList& pointFromDiag,
    const labelListList& pointToFace,
    const label cellID,

    EdgeMap<label>& vertsToFace,
    DynamicList<face>& compactFaces,
    DynamicList<label>& compactCellIDs
) const
{
    // We can form pockets:
    // - 1. triangle on face
    // - 2. multiple triangles on interior (from diag edges)
    // - the edge loop will be pocket since it is only the diag
    //   edges that give it volume?

    // Retriangulate the exterior loops
    const labelListList& edgeLoops = pp.edgeLoops();

    forAll(edgeLoops, loopi)
    {
        const labelList& loop = edgeLoops[loopi];

        if (loop.size() > 2)
        {
            compactFaces.append(face(0));
            face& f = compactFaces.last();

            filterFace
            (
                pp,
                loop,
                pointFromDiag,
                pointToFace,
                f,
                vertsToFace
            );

            compactCellIDs.append(cellID);
        }
    }
}


Foam::MeshedSurface<Foam::face> Foam::isoSurfaceTopo::removeInsidePoints
(
    const bool filterDiag,
    const MeshedSurface<face>& s,
    const boolList& pointFromDiag,
    const labelListList& pointToFace,
    const labelList& start,                 // per cell the starting triangle
    EdgeMap<label>& vertsToFace,
    DynamicList<label>& pointCompactMap,    // per returned point the original
    DynamicList<label>& compactCellIDs      // per returned tri the cellID
) const
{
    const pointField& points = s.points();

    pointCompactMap.clear();
    compactCellIDs.clear();

    DynamicList<face> compactFaces(s.size()/8);

    for (label celli = 0; celli < start.size()-1; celli++)
    {
        // All triangles of the current cell

        label nTris = start[celli+1]-start[celli];

        if (nTris)
        {
            const SubList<face> cellFaces(s, nTris, start[celli]);
            const primitivePatch pp(cellFaces, points);

            triangulateOutside
            (
                filterDiag,
                pp,
                pointFromDiag,
                pointToFace,
                //protectedFace,
                celli,

                vertsToFace,
                compactFaces,
                compactCellIDs
            );
        }
    }


    // Compact out unused points
    // Pick up the used vertices
    labelList oldToCompact(points.size(), -1);
    DynamicField<point> compactPoints(points.size());
    pointCompactMap.clear();

    forAll(compactFaces, i)
    {
        face& f = compactFaces[i];
        forAll(f, fp)
        {
            label pointi = f[fp];
            label compacti = oldToCompact[pointi];
            if (compacti == -1)
            {
                compacti = compactPoints.size();
                oldToCompact[pointi] = compacti;
                compactPoints.append(points[pointi]);
                pointCompactMap.append(pointi);
            }
            f[fp] = compacti;
        }
    }


    MeshedSurface<face> cpSurface;
    cpSurface.reset
    (
        compactPoints.xfer(),
        compactFaces.xfer(),
        xferCopy(s.surfZones())
    );

    return cpSurface;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceTopo::isoSurfaceTopo
(
    const polyMesh& mesh,
    const scalarField& cVals,
    const scalarField& pVals,
    const scalar iso,
    const filterType filter
)
:
    mesh_(mesh),
    cVals_(cVals),
    pVals_(pVals),
    iso_(iso)
{
    if (debug)
    {
        Pout<< "isoSurfaceTopo : iso:" << iso_ << " filter:" << filter << endl;
    }

    tetMatcher tet;

    // Determine if any cut through cell
    List<cellCutType> cellCutTypes;
    const label nCutCells = calcCutTypes(tet, cellCutTypes);

    // Per cell: 5 pyramids cut, each generating 2 triangles
    //  - pointToVerts : from generated iso point to originating mesh verts
    DynamicList<edge> pointToVerts(10*nCutCells);
    //  - pointToFace : from generated iso point to originating mesh face
    DynamicList<labelList> pointToFace(10*nCutCells);
    //  - pointFromDiag : from generated iso point whether is on face diagonal
    DynamicList<bool> pointFromDiag(10*nCutCells);

    // Per cell: number of intersected edges:
    //          - four faces cut so 4 mesh edges + 4 face-diagonal edges
    //          - 4 of the pyramid edges
    EdgeMap<label> vertsToPoint(12*nCutCells);
    DynamicList<label> verts(12*nCutCells);
    // Per cell: 5 pyramids cut (since only one pyramid not cut)
    DynamicList<label> faceLabels(5*nCutCells);
    DynamicList<label> cellLabels(5*nCutCells);


    labelList startTri(mesh_.nCells()+1, 0);

    for (label celli = 0; celli < mesh_.nCells(); celli++)
    {
        startTri[celli] = faceLabels.size();
        if (cellCutTypes[celli] != NOTCUT)
        {
            generateTriPoints
            (
                mesh,
                celli,
                tet.isA(mesh_, celli),

                pointToVerts,
                pointToFace,
                pointFromDiag,

                vertsToPoint,
                verts,
                faceLabels
            );

            for (label i = startTri[celli]; i < faceLabels.size(); i++)
            {
                cellLabels.append(celli);
            }
        }
    }
    startTri[mesh_.nCells()] = faceLabels.size();


    pointToVerts_.transfer(pointToVerts);
    meshCells_.transfer(cellLabels);
    pointToFace_.transfer(pointToFace);

    tmp<pointField> allPoints
    (
        interpolate
        (
            mesh_.cellCentres(),
            mesh_.points()
        )
    );


{
    DynamicList<label> counts;
    forAll(pointToFace_, pointi)
    {
        label nFaces = pointToFace_[pointi].size();

        if (nFaces >= counts.size())
        {
            counts.setSize(nFaces+1, 0);
        }
        counts[nFaces]++;
    }
    DebugVar(counts);
}



    // Assign to MeshedSurface
    faceList allTris(faceLabels.size());
    label verti = 0;
    forAll(allTris, i)
    {
        allTris[i].setSize(3);
        allTris[i][0] = verts[verti++];
        allTris[i][1] = verts[verti++];
        allTris[i][2] = verts[verti++];
    }


    surfZoneList allZones(1);
    allZones[0] = surfZone
    (
        "allFaces",
        allTris.size(),     // size
        0,                  // start
        0                   // index
    );

    MeshedSurface<face>::reset
    (
        allPoints.ref().xfer(),
        allTris.xfer(),
        allZones.xfer()
    );


    // Now:
    // - generated faces and points are assigned to *this
    // - per point we know:
    //  - pointOnDiag: whether it is on a face-diagonal edge
    //  - pointToFace_: from what pyramid (cell+face) it was produced
    //    (note that the pyramid faces are shared between multiple mesh faces)
    //  - pointToVerts_ : originating mesh vertex or cell centre


    if (debug)
    {
        Pout<< "isoSurfaceTopo : generated " << size() << " faces." << endl;
    }


    if (filter != NONE)
    {
        // Triangulate outside (filter edges to cell centres and optionally
        // face diagonals)

        // from generated edge (two vertices) back to originating face
        EdgeMap<label> vertsToFace;
        // Build map from edge vertices to face
        if (filter == DIAGCELL)
        {
            vertsToFace.resize(nCutCells);
            const MeshedSurface<face>& s = *this;
            const pointField& points = s.points();
            for (label celli = 0; celli < startTri.size()-1; celli++)
            {
                // All triangles of the current cell
                label nTris = startTri[celli+1]-startTri[celli];
                if (nTris)
                {
                    const SubList<face> cellFaces(s, nTris, startTri[celli]);
                    const primitivePatch pp(cellFaces, points);

                    // Retriangulate the exterior loops
                    const labelListList& edgeLoops = pp.edgeLoops();

                    forAll(edgeLoops, loopi)
                    {
                        const labelList& loop = edgeLoops[loopi];
                        if (loop.size() > 2)
                        {
                            detectFaceEdges
                            (
                                pp,
                                loop,
                                pointFromDiag,
                                pointToFace,
                                vertsToFace
                            );
                        }
                    }
                }
            }

            // Now vertsToFace holds pairs of surface points that will
            // become an edge on a mesh face. Make sure any 
XXXX Needs indexes of edges on face!!! so we can sync.


        }


    Foam::MeshedSurface<Foam::face> Foam::isoSurfaceTopo::removeInsidePoints
    (
        const bool filterDiag,
        const MeshedSurface<face>& s,
        const boolList& pointFromDiag,
        const labelListList& pointToFace,
        const labelList& start,                 // per cell the starting triangle
        EdgeMap<label>& vertsToFace,



        // back to original point
        DynamicList<label> pointCompactMap(points().size());
        // per returned tri the cellID
        DynamicList<label> compactCellIDs(points().size());
        MeshedSurface<face>::operator=
        (
            removeInsidePoints
            (
                (filter == DIAGCELL ? true : false),
                *this,
                pointFromDiag,
                pointToFace_,
                startTri,
                vertsToFace,
                pointCompactMap,
                compactCellIDs
            )
        );

        pointToVerts_ = UIndirectList<edge>(pointToVerts_, pointCompactMap)();
        pointToFace_ = UIndirectList<labelList>(pointToFace_, pointCompactMap)();
        pointFromDiag = UIndirectList<bool>(pointFromDiag, pointCompactMap)();
        meshCells_.transfer(compactCellIDs);

        if (debug)
        {
            Pout<< "isoSurfaceTopo :"
                << " after removing cell centre and face-diag triangles : "
                << size() << endl;
        }


        if (filter == DIAGCELL)
        {
            // We remove verts on face diagonals. This is in fact just
            // straightening the edges of the face through the cell. This can
            // close off 'pockets' of triangles and create open or
            // multiply-connected triangles

            // Solved by eroding open-edges
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~


            // Mark points on mesh outside. Note that we extend with nCells
            // so we can easily index with pointToVerts_.
            PackedBoolList isBoundaryPoint(mesh.nPoints() + mesh.nCells());
            for
            (
                label facei = mesh.nInternalFaces();
                facei < mesh.nFaces();
                facei++
            )
            {
                isBoundaryPoint.set(mesh.faces()[facei]);
            }


            while (true)
            {
                const labelList& mp = meshPoints();

                PackedBoolList removeFace(this->size());
                label nFaces = 0;
                {
                    const labelListList& edgeFaces =
                        MeshedSurface<face>::edgeFaces();
                    forAll(edgeFaces, edgei)
                    {
                        const labelList& eFaces = edgeFaces[edgei];
                        if (eFaces.size() == 1)
                        {
                            // Open edge. Check that vertices do not originate
                            // from a boundary face
                            const edge& e = edges()[edgei];
                            const edge& verts0 = pointToVerts_[mp[e[0]]];
                            const edge& verts1 = pointToVerts_[mp[e[1]]];
                            if
                            (
                                isBoundaryPoint[verts0[0]]
                             && isBoundaryPoint[verts0[1]]
                             && isBoundaryPoint[verts1[0]]
                             && isBoundaryPoint[verts1[1]]
                            )
                            {
                                // Open edge on boundary face. Keep
                            }
                            else
                            {
                                // Open edge. Mark for erosion
                                if (removeFace.set(eFaces[0]))
                                {
                                    nFaces++;
                                }
                            }
                        }
                    }
                }

                if (debug)
                {
                    Pout<< "isoSurfaceTopo :"
                        << " removing " << nFaces
                        << " faces since on open edges" << endl;
                }

                if (returnReduce(nFaces, sumOp<label>()) == 0)
                {
                    break;
                }

                // Remove the faces
                labelHashSet keepFaces(2*size());
                forAll(removeFace, facei)
                {
                    if (!removeFace[facei])
                    {
                        keepFaces.insert(facei);
                    }
                }

                labelList pointMap;
                labelList faceMap;
                MeshedSurface<face> filteredSurf
                (
                    MeshedSurface<face>::subsetMesh
                    (
                        keepFaces,
                        pointMap,
                        faceMap
                    )
                );
                MeshedSurface<face>::transfer(filteredSurf);

                pointToVerts_ = UIndirectList<edge>(pointToVerts_, pointMap)();
                pointToFace_ = UIndirectList<labelList>(pointToFace_, pointMap)();
                pointFromDiag = UIndirectList<bool>(pointFromDiag, pointMap)();
                meshCells_ = UIndirectList<label>(meshCells_, faceMap)();
            }
        }
    }
}


// ************************************************************************* //
