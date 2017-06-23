/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    levelSetMesh

\*---------------------------------------------------------------------------*/

#include "cut.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "tetIndices.H"
#include "polyMeshTetDecomposition.H"
#include "OBJstream.H"
#include "triSurfaceMesh.H"
#include "mergePoints.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeEdges(OBJstream& str, const FixedList<point, 4>& pts)
{
    str.write(linePointRef(pts[0], pts[1]));
    str.write(linePointRef(pts[1], pts[2]));
    str.write(linePointRef(pts[2], pts[0]));
    str.write(linePointRef(pts[0], pts[3]));
    str.write(linePointRef(pts[1], pts[3]));
    str.write(linePointRef(pts[2], pts[3]));
}


int main(int argc, char *argv[])
{
    argList::validArgs.append("surface");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word surfName = args[1];

    const pointField& cc = mesh.cellCentres();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();


    scalarField levelC;
    scalarField levelP;

    if (false)
    {
        // Distance to plane
        const plane pl
        (
            point(0, 0.0499, 0),
            vector(0, 1, 0)
        );
        Pout<< "refPoint:" << pl.refPoint() << endl;
        Pout<< "normal  :" << pl.normal() << endl;

        // Determine distance to plane
        levelC = ((cc-pl.refPoint()) & pl.normal());
        levelP = ((points-pl.refPoint()) & pl.normal());
    }
    else
    {
        // Distance to surface
        triSurfaceMesh surf
        (
            IOobject
            (
                surfName,
                runTime.constant(),
                "triSurface",
                runTime
            )
        );

        {
            List<pointIndexHit> nearest;
            vectorField normal;
            labelList region;
            surf.searchableSurface::findNearest
            (
                cc,
                scalarField(cc.size(), GREAT),
                nearest,
                normal,
                region
            );

            levelC.setSize(cc.size());
            forAll(levelC, celli)
            {
                vector d(cc[celli]-nearest[celli].hitPoint());
                levelC[celli] = mag(d);

                if ((d&normal[celli]) < 0)
                {
                    levelC[celli] = -levelC[celli];
                }
            }
        }
        {
            List<pointIndexHit> nearest;
            vectorField normal;
            labelList region;
            surf.searchableSurface::findNearest
            (
                points,
                scalarField(points.size(), GREAT),
                nearest,
                normal,
                region
            );

            levelP.setSize(points.size());
            forAll(levelP, pointi)
            {
                vector d(points[pointi]-nearest[pointi].hitPoint());
                levelP[pointi] = mag(d);

                if ((d&normal[pointi]) < 0)
                {
                    levelP[pointi] = -levelP[pointi];
                }
            }
        }
    }


    // Perturb levelC, levelP to be away from 0
    forAll(levelC, i)
    {
        if (mag(levelC[i]) < ROOTSMALL)
        {
            if (levelC[i] >= 0)
            {
                levelC[i] = ROOTSMALL;
            }
            else
            {
                levelC[i] = -ROOTSMALL;
            }
        }
    }
    forAll(levelP, i)
    {
        if (mag(levelP[i]) < ROOTSMALL)
        {
            if (levelP[i] >= 0)
            {
                levelP[i] = ROOTSMALL;
            }
            else
            {
                levelP[i] = -ROOTSMALL;
            }
        }
    }

    DebugVar(levelC);
    DebugVar(levelP);



    typedef DynamicList<FixedList<point, 4>> tetList;
    tetList tets;
    cut::appendOp<tetList> tetAccumulator(tets);

    forAll(cc, cI)
    {
        const List<tetIndices> cellTetIs =
            polyMeshTetDecomposition::cellTetIndices(mesh, cI);

        forAll(cellTetIs, cellTetI)
        {
            const tetIndices& tetIs = cellTetIs[cellTetI];
            const face& f = faces[tetIs.face()];

            const label pI0 = f[tetIs.faceBasePt()];
            const label pIA = f[tetIs.facePtA()];
            const label pIB = f[tetIs.facePtB()];

            const FixedList<point, 4>
                tet =
                {
                    cc[cI],
                    points[pI0],
                    points[pIA],
                    points[pIB]
                };
            const FixedList<scalar, 4>
                level =
                {
                    levelC[cI],
                    levelP[pI0],
                    levelP[pIA],
                    levelP[pIB]
                };

            tetCut(tet, level, tetAccumulator, cut::noOp());
        }
    }

    OBJstream str(runTime.path()/"tets.obj");
    forAll(tets, i)
    {
        writeEdges(str, tets[i]);
    }

    // Construct mesh
    runTime++;

    DynamicList<point> allPoints(4*tets.size());
    forAll(tets, i)
    {
        const FixedList<point, 4>& tet = tets[i];
        allPoints.append(tet[0]);
        allPoints.append(tet[1]);
        allPoints.append(tet[2]);
        allPoints.append(tet[3]);
    }

    labelList oldToNew;
    pointField newPoints;
    mergePoints(allPoints, 100*SMALL, true, oldToNew, newPoints);
    //newPoints = allPoints;
    //oldToNew = identity(allPoints.size());


    cellShapeList cells(tets.size());
    {
        const cellModel& tetShape = *(cellModeller::lookup("tet"));
        labelList tetPoints(4);
        forAll(tets, i)
        {
            tetPoints[0] = oldToNew[4*i];
            tetPoints[1] = oldToNew[4*i+1];
            tetPoints[2] = oldToNew[4*i+2];
            tetPoints[3] = oldToNew[4*i+3];
            cells[i] = cellShape(tetShape, tetPoints, true);
        }
    }

    polyMesh mesh2
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime
        ),
        xferMove(newPoints),
        cells,
        faceListList(0),
        wordList(0),
        wordList(0),
        "defaultFaces",
        polyPatch::typeName,
        wordList(0)
    );

    mesh2.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
