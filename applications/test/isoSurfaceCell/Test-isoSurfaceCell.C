/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Test-isoSurfaceCell

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "fvMesh.H"
#include "isoSurfaceCell.H"
#include "OBJstream.H"
#include "DynamicField.H"
#include "vtkSurfaceWriter.H"
//#include "triSurfaceMesh.H"
#include "smoothTriSurfaceMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "EdgeMap.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// void triangulateOutside
// (
//     const triSurface& s,
//     const label cellID,
// 
//     label start,
//     label size,
// 
//     face& f,
//     DynamicList<face>& triFaces,
//     DynamicList<labelledTri>& compactTris,
//     DynamicList<label>& compactCellIDs
// )
// {
//     const pointField& points = s.points();
// 
//     // All triangles of the current cell
//     SubList<labelledTri> cellTris(s, size, start);
// 
// Pout<< "Triangulating slice size:" << size
//     << " start:" << start << endl;
// 
//     PrimitivePatch<labelledTri, SubList, const pointField&> pp
//     (
//         cellTris,
//         points
//     );
// 
//     // Retriangulate the exterior loops
// 
//     const labelListList& edgeLoops = pp.edgeLoops();
//     const labelList& mp = pp.meshPoints();
// 
//     forAll(edgeLoops, loopi)
//     {
//         const labelList& loop = edgeLoops[loopi];
// 
//         f.setSize(loop.size());
//         forAll(f, i)
//         {
//             f[i] = mp[loop[i]];
//         }
// 
//         triFaces.clear();
//         f.triangles(points, triFaces);
// 
//         forAll(triFaces, i)
//         {
//             const face& f = triFaces[i];
//             compactTris.append(labelledTri(f[0], f[1], f[2], 0));
//             compactCellIDs.append(cellID);
//         }
//     }
// 
//     DebugVar(compactTris);
//     DebugVar(compactCellIDs);
// }
// 
// 
// triSurface removeInsidePoints
// (
//     const triSurface& s,
//     const labelList& cellIDs,
//     const boolList& usesCellCentre,
//     DynamicList<label>& pointCompactMap,    // per returned point the original
//     DynamicList<label>& compactCellIDs      // per returned tri the cellID
// )
// {
//     const pointField& points = s.points();
// 
//     if (cellIDs.size() != s.size() || usesCellCentre.size() != points.size())
//     {
//         FatalErrorInFunction << " Size mismatch" << exit(FatalError);
//     }
// 
//     pointCompactMap.clear();
//     compactCellIDs.clear();
// 
//     DynamicList<face> triFaces;
//     DynamicList<labelledTri> compactTris;
//     face f;
// 
// 
//     label start = 0;
//     forAll(s, trii)
//     {
//         if (trii > 0 && cellIDs[trii] != cellIDs[trii-1])
//         {
//             // All triangles of the current cell
//             triangulateOutside
//             (
//                 s,
//                 cellIDs[trii-1],
// 
//                 start,
//                 trii-start,
// 
//                 f,
//                 triFaces,
//                 compactTris,
//                 compactCellIDs
//             );
// 
//             start = trii;
//         }
//     }
// 
//     // Do final
//     triangulateOutside
//     (
//         s,
//         cellIDs[cellIDs.size()-1],
// 
//         start,
//         cellIDs.size()-start,
// 
//         f,
//         triFaces,
//         compactTris,
//         compactCellIDs
//     );
// 
// 
// 
// 
//     // Compact out unused points
//     // Pick up the used vertices
//     labelList oldToCompact(points.size(), -1);
//     DynamicField<point> compactPoints(points.size());
//     pointCompactMap.clear();
// 
//     forAll(compactTris, i)
//     {
//         labelledTri& f = compactTris[i];
//         forAll(f, fp)
//         {
//             label pointi = f[fp];
//             label compacti = oldToCompact[pointi];
//             if (compacti == -1)
//             {
//                 compacti = compactPoints.size();
//                 oldToCompact[pointi] = compacti;
//                 compactPoints.append(points[pointi]);
//                 pointCompactMap.append(pointi);
//             }
//             f[fp] = compacti;
//         }
//     }
// 
//     return triSurface
//     (
//         compactTris.xfer(),
//         geometricSurfacePatchList(0),
//         compactPoints.xfer()
//     );
// }
// 
// 
// void filterDiag
// (
//     const isoSurfaceCell& iso,
//     const label start,
//     const label size,
//     DynamicList<face>& faces,
//     DynamicField<scalar>& faceMeshCells
// )
// {
//     const labelList& meshCells = iso.meshCells();
// 
//     // All triangles of the current cell
//     SubList<labelledTri> cellTris(iso, size, start);
// 
//     PrimitivePatch<labelledTri, SubList, const pointField&> pp
//     (
//         cellTris,
//         iso.points()
//     );
// 
//     const labelListList& edgeLoops = pp.edgeLoops();
//     const labelList& mp = pp.meshPoints();
// 
//     forAll(edgeLoops, loopi)
//     {
//         const labelList& loop = edgeLoops[loopi];
// 
//         bool filtered = false;
//         if (loop.size() > 2)
//         {
//             face f(loop.size());
//             label fpi = 0;
//             forAll(f, i)
//             {
//                 if (!iso.isOnDiag()[mp[loop[i]]])
//                 {
//                     f[fpi++] = mp[loop[i]];
//                 }
//             }
//             if (fpi > 2)
//             {
//                 f.setSize(fpi);
//                 faces.append(f);
//                 faceMeshCells.append(scalar(meshCells[start]));
//                 filtered = true;
//             }
//             else
//             {
//                 Pout<< "Discarding filtered face:"
//                     << SubList<label>(f, fpi)
//                     << " for cell:" << SubList<label>(meshCells, size, start)
//                     << " with triangles:" << pp
//                     << " with trianglescc:" << pp.faceCentres()
//                     << endl;
//             }   
//         }
// 
//         if (!filtered)
//         {
//             faces.append(face(UIndirectList<label>(mp, loop)()));
//             faceMeshCells.append(scalar(meshCells[start]));
// 
// Pout<< "For cell:" << SubList<label>(meshCells, size, start)
//     << " have triangles:" << pp
//     << " have trianglescc:" << pp.faceCentres()
//     << " have face:" << faces.last() << endl;
//         }
//     }
// }


label generatePoint
(
    const scalar isoValue,
    const scalar v0,
    const point& pt0,
    const label i0,                     // index

    const scalar v1,
    const point& pt1,
    const label i1,                     // index

    DynamicList<point>& points,
    EdgeMap<label>& endToPoint
)
{
    scalar s = (isoValue-v0)/(v1-v0);

    EdgeMap<label>::const_iterator edgeFnd = endToPoint.find(edge(i0, i1));
    if (edgeFnd != endToPoint.end())
    {
        Pout<< "For edge between " << pt0
            << " and " << pt1 << " found existing point " << edgeFnd()
            << " at location:" << points[edgeFnd()] << endl;
        return edgeFnd();
    }
    else
    {
        // Generate new point
        label pointi = points.size();
        points.append(s*pt1+(1.0-s)*pt0);
        endToPoint.insert(edge(i0, i1), pointi);
        return pointi;
    }
}
void generateEdgeCut
(
    const scalar isoValue,
    const scalar v0,
    const point& pt0,
    const label i0,                     // index

    const scalar v1,
    const point& pt1,
    const label i1,                     // index

    DynamicList<point>& points,
    EdgeMap<label>& endToPoint,
    DynamicList<label>& verts
)
{
    if
    (
        (v0 < isoValue && v1 > isoValue)
     || (v0 > isoValue && v1 < isoValue)
    )
    {
        verts.append
        (
            generatePoint
            (
                isoValue, v0, pt0, i0, v1, pt1, i1, 
                points,
                endToPoint
            )
        );
    }
}


void generateEdgeCuts
(
    const scalar isoValue,
    const scalar v0,
    const point& pt0,
    const label i0,                     // index

    const scalar v1,
    const point& pt1,
    const label i1,                     // index

    const scalar v2,
    const point& pt2,
    const label i2,                     // index

    const scalar v3,
    const point& pt3,
    const label i3,

    DynamicList<point>& points,
    EdgeMap<label>& table,
    DynamicList<label>& verts
)
{
    generateEdgeCut(isoValue, v0, pt0, i0, v1, pt1, i1, points, table, verts);
    generateEdgeCut(isoValue, v0, pt0, i0, v2, pt2, i2, points, table, verts);
    generateEdgeCut(isoValue, v0, pt0, i0, v3, pt3, i3, points, table, verts);
    generateEdgeCut(isoValue, v1, pt1, i1, v2, pt2, i2, points, table, verts);
    generateEdgeCut(isoValue, v2, pt2, i2, v3, pt3, i3, points, table, verts);
}
void generateTriPoints
(
    const scalar isoValue,
    const FixedList<scalar, 4>& s,
    const FixedList<point, 4>& p,
    const FixedList<label, 4>& pIndex,

    DynamicList<point>& points,
    EdgeMap<label>& endToPoint,
    DynamicList<label>& verts       // every three verts is new triangle
)
{
    int triIndex = 0;
    if (s[0] < isoValue)
    {
        triIndex |= 1;
    }
    if (s[1] < isoValue)
    {
        triIndex |= 2;
    }
    if (s[2] < isoValue)
    {
        triIndex |= 4;
    }
    if (s[3] < isoValue)
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
                    isoValue,
                    s[0], p[0], pIndex[0], s[1], p[1], pIndex[1],
                    points,
                    endToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[0],p[0],pIndex[0],s[2],p[2],pIndex[2],
                    points,
                    endToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[0],p[0],pIndex[0],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
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
                    isoValue,
                    s[1],p[1],pIndex[1],s[0],p[0],pIndex[0],
                    points,
                    endToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[1],p[1],pIndex[1],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[1],p[1],pIndex[1],s[2],p[2],pIndex[2],
                    points,
                    endToPoint
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
                    isoValue,
                    s[0],p[0],pIndex[0],s[2],p[2],pIndex[2],
                    points,
                    endToPoint
                )
            );
            label p1p3
            (
                generatePoint
                (
                    isoValue,
                    s[1],p[1],pIndex[1],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
                )
            );

            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[0],p[0],pIndex[0],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
                )
            );
            verts.append(p1p3);
            verts.append(p0p2);
            verts.append(p1p3);
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[1],p[1],pIndex[1],s[2],p[2],pIndex[2],
                    points,
                    endToPoint
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
                    isoValue,
                    s[2],p[2],pIndex[2],s[0],p[0],pIndex[0],
                    points,
                    endToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[2],p[2],pIndex[2],s[1],p[1],pIndex[1],
                    points,
                    endToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[2],p[2],pIndex[2],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
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
                    isoValue,
                    s[0],p[0],pIndex[0],s[1],p[1],pIndex[1],
                    points,
                    endToPoint
                )
            );
            label p2p3
            (
                generatePoint
                (
                    isoValue,
                    s[2],p[2],pIndex[2],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
                )
            );

            verts.append(p0p1);
            verts.append(p2p3);
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[0],p[0],pIndex[0],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
                )
            );
            verts.append(p0p1);
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[1],p[1],pIndex[1],s[2],p[2],pIndex[2],
                    points,
                    endToPoint
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
                    isoValue,
                    s[0],p[0],pIndex[0],s[1],p[1],pIndex[1],
                    points,
                    endToPoint
                )
            );
            label p2p3
            (
                generatePoint
                (
                    isoValue,
                    s[2],p[2],pIndex[2],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
                )
            );

            verts.append(p0p1);
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[1],p[1],pIndex[1],s[3],p[3],pIndex[3],
                    points,
                    endToPoint
                )
            );
            verts.append(p2p3);
            verts.append(p0p1);
            verts.append(p2p3);
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[0],p[0],pIndex[0],s[2],p[2],pIndex[2],
                    points,
                    endToPoint
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
                    isoValue,
                    s[3],p[3],pIndex[3],s[0],p[0],pIndex[0],
                    points,
                    endToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[3],p[3],pIndex[3],s[2],p[2],pIndex[2],
                    points,
                    endToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    isoValue,
                    s[3],p[3],pIndex[3],s[1],p[1],pIndex[1],
                    points,
                    endToPoint
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



// isosurface algorithm
void generateEdgeCuts
(
    const polyMesh& mesh,
    const scalar isoValue,
    const scalarField& cellValues,
    const scalarField& pointValues,
    const label celli,

    DynamicList<point>& points,
    EdgeMap<label>& endToPoint,
    DynamicList<label>& verts
)
{
    const cell& cFaces = mesh.cells()[celli];
    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];
        const face& f = mesh.faces()[facei];

        label fp0 = mesh.tetBasePtIs()[facei];

        // Skip undefined tets
        if (fp0 < 0)
        {
            fp0 = 0;
        }

        label fp = f.fcIndex(fp0);
        for (label i = 2; i < f.size(); i++)
        {
            label nextFp = f.fcIndex(fp);

            label p0 = f[fp0];
            label p1 = f[fp];
            label p2 = f[nextFp];
            if (mesh.faceOwner()[facei] == celli)
            {
                Swap(p1, p2);
            }

            generateTriPoints
            (
                isoValue,
                FixedList<scalar, 4>
                ({
                    pointValues[p0],
                    pointValues[p1],
                    pointValues[p2],
                    cellValues[celli]
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

                points,
                endToPoint,
                verts       // every three verts is new triangle
            );

            // generateEdgeCuts
            // (
            //     isoValue,
            //     p0,                     // index
            //     pointValues[p0],
            //     mesh.points()[p0],
            //     p1,                     // index
            //     pointValues[p1],
            //     mesh.points()[p1],
            //     p2,                     // index
            //     pointValues[p2],
            //     mesh.points()[p2],
            // 
            //     mesh.nPoints()+celli,   // index
            //     cellValues[celli],
            //     mesh.cellCentres()[celli],
            // 
            //     points,
            //     endToPoint,
            //     verts
            // );

            fp = nextFp;
        }
    }
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    Random rndGen(0);

    scalarField cellValues;
    scalarField pointValues;
    scalar isoValue;

    // Random or from positions
    if (true)
    {
        cellValues = mesh.cellCentres().component(vector::Y);
        pointValues = mesh.points().component(vector::Y);
        //const scalar minPoints(min(pointValues));
        //const scalar maxPoints(max(pointValues));
        //forAll(pointValues, i)
        //{
        //    pointValues[i] = minPoints+(maxPoints-minPoints)*rndGen.scalar01();
        //}
        isoValue = 0.51*(average(cellValues)+average(pointValues));

        DynamicList<point> points;
        EdgeMap<label> endToPoint;
        DynamicList<label> verts;
        generateEdgeCuts
        (
            mesh,
            isoValue,
            cellValues,
            pointValues,
            0,

            points,
            endToPoint,
            verts
        );

        List<labelledTri> tris(verts.size()/3);
        label verti = 0;
        forAll(tris, i)
        {
            label v0 = verts[verti++];
            label v1 = verts[verti++];
            label v2 = verts[verti++];
            tris[i] = labelledTri(v0, v1, v2, 0);
        }
        triSurface s(tris, pointField(points));
        s.write("simple.obj");
    }
    if (true)
    {
        const smoothTriSurfaceMesh searchSurf
        (
            IOobject
            (
                "object.stl",
                runTime.constant(),
                "triSurface",
                runTime
            ),
            180
        );
        List<pointIndexHit> cellNearest;
        vectorField cellNormal;
        {
            const pointField& cc = mesh.cellCentres();

            searchSurf.findNearest
            (
                cc,
                scalarField(cc.size(), great),
                cellNearest
            );
            searchSurf.getNormal(cellNearest, cellNormal);

            cellValues.setSize(cc.size());
            forAll(cc, celli)
            {
                const vector d(cc[celli]-cellNearest[celli].hitPoint());
                cellValues[celli] = sign(d&cellNormal[celli])*mag(d);
            }
        }
        List<pointIndexHit> pointNearest;
        vectorField pointNormal;
        {
            const pointField& points = mesh.points();

            searchSurf.findNearest
            (
                points,
                scalarField(points.size(), great),
                pointNearest
            );
            searchSurf.getNormal(pointNearest, pointNormal);

            pointValues.setSize(points.size());
            forAll(points, pointi)
            {
                const vector d(points[pointi]-pointNearest[pointi].hitPoint());
                pointValues[pointi] = sign(d&pointNormal[pointi])*mag(d);
            }

            isoValue = 0.0;
        }


        // Check inconsistency: change of sign inside cell & cell further
        // away
        OBJstream farCells(runTime.path()/"straddlingCells.obj");
        forAll(cellNearest, celli)
        {
            const point& cc = mesh.cellCentres()[celli];
            const point& near = cellNearest[celli].hitPoint();
            if (mag(near-cc) > 0.1)
            {
                const labelList& cPoints = mesh.cellPoints()[celli];
                bool straddle = false;
                forAll(cPoints, i)
                {
                    label pointi = cPoints[i];
                    if (sign(cellValues[celli]) != sign(pointValues[pointi]))
                    {
                        straddle = true;
                        break;
                    }
                }
                if (straddle)
                {
                    Pout<< "DEtected cell " << cc
                        << " with multiple signs." << endl;

                    const cell& cFaces = mesh.cells()[celli];
                    forAll(cFaces, cFacei)
                    {
                        const face& f = mesh.faces()[cFaces[cFacei]];
                        farCells.write(f, mesh.points());
                    }

                    farCells.write(linePointRef(cc, near));
                    forAll(cPoints, i)
                    {
                        label pointi = cPoints[i];
                        const point& pt = mesh.points()[pointi];
                        const point& nearPt = pointNearest[pointi].hitPoint();
                        farCells.write(linePointRef(pt, nearPt));
                    }
                }
            }
        }
        Pout<< "Written " << farCells.nVertices() << " to " << farCells.name()
            << endl;
    }


    // Optimisation 1: isoSurfaceCell can remove the cell centre

    //- From point field and interpolated cell.
    scalarField cellAvg(mesh.nCells(), scalar(0));
    {
        labelField nPointCells(mesh.nCells(), 0);
        {
            for (label pointi = 0; pointi < mesh.nPoints(); pointi++)
            {
                const labelList& pCells = mesh.pointCells(pointi);

                forAll(pCells, i)
                {
                    label celli = pCells[i];
                    cellAvg[celli] += pointValues[pointi];
                    nPointCells[celli]++;
                }
            }
        }
        forAll(cellAvg, celli)
        {
            cellAvg[celli] /= nPointCells[celli];
        }
    }


    {
        volScalarField cDist
        (
            IOobject
            (
                "cellDistance",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimLength, 0)
        );
        cDist.primitiveFieldRef()= cellValues;
        cDist.correctBoundaryConditions();
        cDist.write();

        pointScalarField pDist
        (
            IOobject
            (
                "pointDistance",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(mesh),
            dimensionedScalar("zero", dimLength, 0)
        );
        pDist.primitiveFieldRef() = pointValues;
        pDist.write();
    }

    isoSurfaceCell iso
    (
        mesh,
        cellValues,     //cellAvg
        pointValues,
        isoValue,
        true       // regularise = remove cell centre
    );

    Pout<< "iso:" << iso.size() << endl;

    iso.write("iso.obj");

    {
        vtkSurfaceWriter vtk;
        const labelList& meshCells = iso.meshCells();
        faceList localFaces(meshCells.size());
        scalarField scalarMeshCells(meshCells.size());
        forAll(meshCells, i)
        {
            localFaces[i] = face(iso.localFaces()[i]);
            scalarMeshCells[i] = 1.0*meshCells[i];
        }
        vtk.write
        (
            runTime.path(),
            "iso_with_meshCells",
            iso.localPoints(),  //iso.points(),
            localFaces,   //faces,
            "meshCells",
            scalarMeshCells,
            false
        );
    }

// 
//     {
//         OBJstream str("isOnDiag.obj");
//         forAll(iso.isOnDiag(), pointi)
//         {
//             if (iso.isOnDiag()[pointi])
//             {
//                 str.write(iso.points()[pointi]);
//             }
//         }
//     }
// 
// 
//     // Optimisation 2: keep the outside loop only. Make a single face
//     // of all of the interior
// 
// 
//     // Collect triangles per face and filter according to cells
//     const labelList& meshCells = iso.meshCells();
// 
// 
//     DynamicList<face> faces;
//     DynamicField<scalar> faceMeshCells;
//     {
//         label start = 0;
//         forAll(meshCells, trii)
//         {
//             if (meshCells[trii] != meshCells[start])
//             {
//                 filterDiag
//                 (
//                     iso,
//                     start,
//                     trii-start,
//                     faces,
//                     faceMeshCells
//                 );
// 
//                 start = trii;
//             }
//         }
//         // Do the last ones
//         if (meshCells.size())
//         {
//             filterDiag
//             (
//                 iso,
//                 start,
//                 meshCells.size()-start,
//                 faces,
//                 faceMeshCells
//             );
//         }
//     }
// 
//     primitiveFacePatch compact(faces, iso.points());
// 
//     vtkSurfaceWriter vtk;
//     vtk.write
//     (
//         runTime.path(),
//         "mySurface",
//         compact.localPoints(),  //iso.points(),
//         compact.localFaces(),   //faces,
//         "meshCells",
//         faceMeshCells,
//         false
//     );

//    OBJstream str("ccPoints.obj");
//    forAll(iso.usesCellCentre(), pointi)
//    {
//        if (iso.usesCellCentre()[pointi])
//        {
//            str.write(iso.points()[pointi]);
//        }
//    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
