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
#include "triSurfaceMesh.H"
#include "volFields.H"
#include "pointFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triangulateOutside
(
    const triSurface& s,
    const label cellID,

    label start,
    label size,

    face& f,
    DynamicList<face>& triFaces,
    DynamicList<labelledTri>& compactTris,
    DynamicList<label>& compactCellIDs
)
{
    const pointField& points = s.points();

    // All triangles of the current cell
    SubList<labelledTri> cellTris(s, size, start);

Pout<< "Triangulating slice size:" << size
    << " start:" << start << endl;

    PrimitivePatch<labelledTri, SubList, const pointField&> pp
    (
        cellTris,
        points
    );

    // Retriangulate the exterior loops

    const labelListList& edgeLoops = pp.edgeLoops();
    const labelList& mp = pp.meshPoints();

    forAll(edgeLoops, loopi)
    {
        const labelList& loop = edgeLoops[loopi];

        f.setSize(loop.size());
        forAll(f, i)
        {
            f[i] = mp[loop[i]];
        }

        triFaces.clear();
        f.triangles(points, triFaces);

        forAll(triFaces, i)
        {
            const face& f = triFaces[i];
            compactTris.append(labelledTri(f[0], f[1], f[2], 0));
            compactCellIDs.append(cellID);
        }
    }

    DebugVar(compactTris);
    DebugVar(compactCellIDs);
}


triSurface removeInsidePoints
(
    const triSurface& s,
    const labelList& cellIDs,
    const boolList& usesCellCentre,
    DynamicList<label>& pointCompactMap,    // per returned point the original
    DynamicList<label>& compactCellIDs      // per returned tri the cellID
)
{
    const pointField& points = s.points();

    if (cellIDs.size() != s.size() || usesCellCentre.size() != points.size())
    {
        FatalErrorInFunction << " Size mismatch" << exit(FatalError);
    }

    pointCompactMap.clear();
    compactCellIDs.clear();

    DynamicList<face> triFaces;
    DynamicList<labelledTri> compactTris;
    face f;


    label start = 0;
    forAll(s, trii)
    {
        if (trii > 0 && cellIDs[trii] != cellIDs[trii-1])
        {
            // All triangles of the current cell
            triangulateOutside
            (
                s,
                cellIDs[trii-1],

                start,
                trii-start,

                f,
                triFaces,
                compactTris,
                compactCellIDs
            );

            start = trii;
        }
    }

    // Do final
    triangulateOutside
    (
        s,
        cellIDs[cellIDs.size()-1],

        start,
        cellIDs.size()-start,

        f,
        triFaces,
        compactTris,
        compactCellIDs
    );




    // Compact out unused points
    // Pick up the used vertices
    labelList oldToCompact(points.size(), -1);
    DynamicField<point> compactPoints(points.size());
    pointCompactMap.clear();

    forAll(compactTris, i)
    {
        labelledTri& f = compactTris[i];
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

    return triSurface
    (
        compactTris.xfer(),
        geometricSurfacePatchList(0),
        compactPoints.xfer()
    );
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
    if (false)
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
    }
    if (true)
    {
        const triSurfaceMesh searchSurf
        (
            IOobject
            (
                "object.stl",
                runTime.constant(),
                "triSurface",
                runTime
            )
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
    }

//    {
//        OBJstream farCells(runTime.path()/"straddlingCells.obj");
//
//        const labelListList& cellPoints = mesh.cellPoints();
//        forAll(cellPoints, celli)
//        {
//            bool cSign = sign(cellValues[celli]);
//        }
//    }



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
        false       // regularise = remove cell centre
    );

    Pout<< "iso:" << iso << endl;

    iso.write("iso.obj");

//    {
//        OBJstream str("isOnDiag.obj");
//        forAll(iso.isOnDiag(), pointi)
//        {
//            if (iso.isOnDiag()[pointi])
//            {
//                str.write(iso.points()[pointi]);
//            }
//        }
//    }


    // Optimisation 2: keep the outside loop only. Make a single face
    // of all of the interior


    // Collect triangles per face and filter according to cells
    const labelList& meshCells = iso.meshCells();


    DynamicList<face> faces;
    DynamicField<scalar> faceMeshCells;
    {
        label start = 0;
        forAll(meshCells, trii)
        {
            if (meshCells[trii] != meshCells[start])
            {
                // All triangles of the current cell
                SubList<labelledTri> cellTris(iso, trii-start, start);

                PrimitivePatch<labelledTri, SubList, const pointField&> pp
                (
                    cellTris,
                    iso.points()
                );

                {
                    const labelListList& edgeLoops = pp.edgeLoops();
                    const labelList& mp = pp.meshPoints();

                    forAll(edgeLoops, loopi)
                    {
                        const labelList& loop = edgeLoops[loopi];

                        if (loop.size() > 2)
                        {
                            face f(loop.size());
                            label fpi = 0;
                            forAll(f, i)
                            {
                                if (!iso.isOnDiag()[mp[loop[i]]])
                                {
                                    f[fpi++] = mp[loop[i]];
                                }
                            }
                            if (fpi > 2)
                            {
                                f.setSize(fpi);
                                faces.append(f);
                                faceMeshCells.append(scalar(meshCells[start]));
                            }
                        }
                    }
                }


                start = trii;
            }
        }
        // Do the last ones
        if (meshCells.size())
        {
            // All triangles of the current cell
            SubList<labelledTri> cellTris(iso, meshCells.size()-start, start);

            PrimitivePatch<labelledTri, SubList, const pointField&> pp
            (
                cellTris,
                iso.points()
            );

            {
                const labelListList& edgeLoops = pp.edgeLoops();
                const labelList& mp = pp.meshPoints();

                forAll(edgeLoops, loopi)
                {
                    const labelList& loop = edgeLoops[loopi];

                    if (loop.size() > 2)
                    {
                        face f(loop.size());
                        label fpi = 0;
                        forAll(f, i)
                        {
                            if (!iso.isOnDiag()[mp[loop[i]]])
                            {
                                f[fpi++] = mp[loop[i]];
                            }
                        }
                        if (fpi > 2)
                        {
                            f.setSize(fpi);
                            faces.append(f);
                            faceMeshCells.append(scalar(meshCells[start]));
                        }
                    }
                }
            }
        }
    }

    primitiveFacePatch compact(faces, iso.points());

    vtkSurfaceWriter vtk;
    vtk.write
    (
        runTime.path(),
        "mySurface",
        compact.localPoints(),  //iso.points(),
        compact.localFaces(),   //faces,
        "meshCells",
        faceMeshCells,
        false
    );

//    OBJstream str("ccPoints.obj");
//    forAll(iso.usesCellCentre(), pointi)
//    {
//        if (iso.usesCellCentre()[pointi])
//        {
//            str.write(iso.points()[pointi]);
//        }
//    }


//    {
//        const labelListList& edgeLoops = iso.edgeLoops();
//        const labelList& mp = iso.meshPoints();
//
//        DynamicList<face> triFaces;
//        face f;
//        forAll(edgeLoops, loopi)
//        {
//            const labelList& loop = edgeLoops[loopi];
//
//            f.setSize(loop.size());
//            forAll(f, i)
//            {
//                f[i] = mp[loop[i]];
//            }
//            f.triangles(iso.points(), triFaces);
//        }
//
//
//        // Pick up the used vertices
//        labelList oldToCompact(mp.size(), -1);
//        DynamicField<point> compactPoints(mp.size());
//
//        forAll(triFaces, i)
//        {
//            face& f = triFaces[i];
//            forAll(f, fp)
//            {
//                label pointi = f[fp];
//                label compacti = oldToCompact[pointi];
//                if (compacti == -1)
//                {
//                    compacti = compactPoints.size();
//                    oldToCompact[pointi] = compacti;
//                    compactPoints.append(iso.points()[pointi]);
//                }
//                f[fp] = compacti;
//            }
//        }
//
//
//        OBJstream str("loops.obj");
//        str.write(triFaces.shrink(), compactPoints, false);
//    }
//    {
//        DynamicList<label> pointCompactMap;
//        DynamicList<label> compactCellIDs;
//        triSurface newSurf
//        (
//            removeInsidePoints
//            (
//                iso,
//                iso.meshCells(),
//                iso.usesCellCentre(),
//                pointCompactMap,    // per returned point the original
//                compactCellIDs      // per returned tri the cellID
//            )
//        );
//
//
//        OBJstream str("compact.obj");
//        newSurf.write("compact.obj");
//    }


//    point avgVert(point::zero);
//    label nVert = 0;
//
//    forAll(iso.usesCellCentre(), pointi)
//    {
//        if (!iso.usesCellCentre()[pointi])
//        {
//            avgVert += iso.points()[pointi];
//            nVert++;
//        }
//    }
//
//    avgVert /= nVert;
//
//    label nCc = iso.points().size()-nVert;
//    if (nCc > 0)
//    {
//        pointField newPoints(iso.points());
//        forAll(iso.usesCellCentre(), pointi)
//        {
//            if (iso.usesCellCentre()[pointi])
//            {
//                newPoints[pointi] = avgVert;
//            }
//        }
//        iso.movePoints(newPoints);
//        iso.write("newIso.obj");
//    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
