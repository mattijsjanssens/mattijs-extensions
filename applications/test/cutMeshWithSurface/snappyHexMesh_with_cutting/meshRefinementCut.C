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

#include "meshRefinement.H"
#include "refinementSurfaces.H"
#include "DynamicField.H"
#include "cellCuts.H"
#include "meshCutter.H"
#include "polyTopoChange.H"
#include "volumeType.H"
#include "searchableSurfaces.H"
#include "OBJstream.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshRefinement::snapToSurface
(
    labelList& pointSurfaceRegion,
    labelList& edgeSurfaceRegion,
    scalarField& edgeWeight
)
{
    const edgeList& edges = mesh_.edges();
    const pointField& points = mesh_.points();

    pointSurfaceRegion.setSize(points.size());
    pointSurfaceRegion = -1;

    edgeSurfaceRegion.setSize(edges.size());
    edgeSurfaceRegion = -1;

    edgeWeight.setSize(edges.size());
    edgeWeight = -GREAT;


    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    //vectorField normal1;

    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    //vectorField normal2;
    {
        vectorField start(edges.size());
        vectorField end(edges.size());
        forAll(edges, edgei)
        {
            const edge& e = edges[edgei];
            start[edgei] = points[e[0]];
            end[edgei] = points[e[1]];
        }
        
        surfaces_.findNearestIntersection
        (
            //labelList(1, 0),    //identity(surfaces_.surfaces().size()),
            identity(surfaces_.surfaces().size()),
            start,
            end,

            surface1,
            hit1,
            region1,
            //normal1,

            surface2,
            hit2,
            region2
            //normal2
        );
    }

    // Adjust location
    // ~~~~~~~~~~~~~~~

    pointField newPoints(points);
    label nAdjusted = 0;

    const labelListList& pointEdges = mesh_.pointEdges();
    forAll(pointEdges, pointi)
    {
        const point& pt = points[pointi];
        const labelList& pEdges = pointEdges[pointi];

        // Get the nearest intersection
        label minEdgei = -1;
        scalar minFraction = 0.5;   // Harpoon 0.25; // Samm?
        forAll(pEdges, pEdgei)
        {
            label edgei = pEdges[pEdgei];
            if (hit1[edgei].hit())
            {
                const point& hitPt = hit1[edgei].hitPoint();

                const edge& e = edges[edgei];
                label otherPointi = e.otherVertex(pointi);
                const point& otherPt = points[otherPointi];

                vector eVec(otherPt-pt);
                scalar f = eVec&(hitPt-pt)/magSqr(eVec);

                if (f < minFraction)
                {
                    minEdgei = edgei;
                    minFraction = f;
                }
            }
        }
        if (minEdgei != -1 && minFraction >= 0.01)
        {
            // Move point to intersection with minEdgei
            if (pointSurfaceRegion[pointi] == -1)
            {
                pointSurfaceRegion[pointi] = surfaces_.globalRegion
                (
                    surface1[minEdgei],
                    region1[minEdgei]
                );
                newPoints[pointi] = hit1[minEdgei].hitPoint();
                nAdjusted++;
            }
        }
    }

    reduce(nAdjusted, sumOp<label>());

    Info<< "Snapped " << nAdjusted << " points" << endl;

    if (nAdjusted > 0)
    {
        mesh_.movePoints(newPoints);

        // Delete mesh volumes and meshPhi. No other way to do this?
        //mesh_.clearOut();

        // Update the intersections on all the moved points
        {
            const labelListList& pointFaces = mesh_.pointFaces();

            PackedBoolList isChangedFace(mesh_.nFaces());
            forAll(pointSurfaceRegion, pointi)
            {
                if (pointSurfaceRegion[pointi] != -1)
                {
                    isChangedFace.set(pointFaces[pointi]);
                }
            }
            DynamicList<label> changedFaces(isChangedFace.count());
            forAll(isChangedFace, facei)
            {
                if (isChangedFace[facei])
                {
                    changedFaces.append(facei);
                }
            }
            updateIntersections(changedFaces);
        }


//DebugVar(hit1.size());
//DebugVar(hit2.size());
//DebugVar(edges.size());

        // Filter out edge cuts if endpoints have been cut
        forAll(hit1, edgei)
        {
            if (hit1[edgei].hit())
            {
                const edge& e = edges[edgei];

                if
                (
                    pointSurfaceRegion[e[0]] == -1
                 && pointSurfaceRegion[e[1]] == -1
                )
                {
                    const point& pt = hit1[edgei].hitPoint();

                    vector eVec(e.vec(points));
                    scalar f = eVec&(pt-points[e[0]])/magSqr(eVec);

                    if (f < 0.01)
                    {
                        pointSurfaceRegion[e[0]] = surfaces_.globalRegion
                        (
                            surface1[edgei],
                            region1[edgei]
                        );
                    }
                    else if (f > 0.99)
                    {
                        pointSurfaceRegion[e[1]] = surfaces_.globalRegion
                        (
                            surface1[edgei],
                            region1[edgei]
                        );
                    }
                }
            }
        }

        forAll(hit1, edgei)
        {
            if (hit1[edgei].hit())
            {
                const edge& e = edges[edgei];

                if
                (
                    pointSurfaceRegion[e[0]] == -1
                 && pointSurfaceRegion[e[1]] == -1
                )
                {
                    const point& pt = hit1[edgei].hitPoint();
                    vector eVec(e.vec(points));
                    scalar f = eVec&(pt-points[e[0]])/magSqr(eVec);

                    edgeSurfaceRegion[edgei] = surfaces_.globalRegion
                    (
                        surface1[edgei],
                        region1[edgei]
                    );
                    edgeWeight[edgei] = f;
                }
            }
        }



        if (debug&MESH)
        {
            const_cast<Time&>(mesh_.time())++;
            Pout<< "Writing snapped mesh to time " << timeName()
                << endl;

            // Prevent meshPhi from being written.
            mesh_.clearOut();

            write
            (
                debugType(debug),
                writeType(writeLevel() | WRITEMESH),
                mesh_.time().path()/"snapped"
            );
        }
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::cutWithSurface
(
    const labelList& pointSurfaceRegion,
    const labelList& edgeSurfaceRegion,
    const scalarField& edgeWeight
)
{
    DynamicList<label> cutVerts;
    DynamicList<label> cutEdges;
    DynamicField<scalar> cutEdgeWeights;
    {
        label n = 0;
        forAll(pointSurfaceRegion, pointi)
        {
            if (pointSurfaceRegion[pointi] != -1)
            {
                n++;
            }
        }
        cutVerts.setCapacity(n);
        forAll(pointSurfaceRegion, pointi)
        {
            if (pointSurfaceRegion[pointi] != -1)
            {
                cutVerts.append(pointi);
            }
        }

        n = 0;
        forAll(edgeSurfaceRegion, edgei)
        {
            if (edgeSurfaceRegion[edgei] != -1)
            {
                n++;
            }
        }
        cutEdges.setCapacity(n);
        cutEdgeWeights.setCapacity(n);
        forAll(edgeSurfaceRegion, edgei)
        {
            if (edgeSurfaceRegion[edgei] != -1)
            {
                cutEdges.append(edgei);
                cutEdgeWeights.append(edgeWeight[edgei]);
            }
        }
    }

    // Write some obj file
    {
        const edgeList& edges = mesh_.edges();
        const pointField& points = mesh_.points();

        OBJstream str(mesh_.time().path()/"allCuts.obj");
        Pout<< "Writing " << returnReduce(cutVerts.size(), sumOp<label>())
            << " cut vertices and "
            << returnReduce(cutEdges.size(), sumOp<label>())
            << " cut edges to " << str.name() << endl;
        forAll(cutVerts, i)
        {
            str.write(points[cutVerts[i]]);
        }
        forAll(cutEdges, i)
        {
            const edge& e = edges[cutEdges[i]];
            const scalar f = cutEdgeWeights[i];
            str.write((1.0-f)*points[e[0]]+f*points[e[1]]);
        }
    }



    //- Construct from pattern of cuts. Detect cells to cut.
    cellCuts cuts(mesh_, cutVerts, cutEdges, cutEdgeWeights);

    DebugVar(cuts.nLoops());

    // See if we need to flip anything
    {
        const pointField& points = mesh_.points();
        const labelListList& cellLoops = cuts.cellLoops();
        const labelListList& cellAnchorPoints = cuts.cellAnchorPoints();


        DynamicList<label> testCells(cuts.nLoops());
        DynamicField<point> testPoints(cuts.nLoops());

        forAll(cellLoops, celli)
        {
            if (cellLoops[celli].size())
            {
                const labelList& anchors = cellAnchorPoints[celli];
                testCells.append(celli);
                testPoints.append(points[anchors[0]]);
            }
        }


        labelList insideSurfaces;

        // Replacement for refinementSurfaces::findInside
        {
            insideSurfaces.setSize(testPoints.size());
            insideSurfaces = -1;

            forAll(surfaces_.surfaces(), surfi)
            {
                label geomi = surfaces_.surfaces()[surfi];
                //label geomi = 0;
                const searchableSurface& geom = surfaces_.geometry()[geomi];

                if (geom.hasVolumeType())
                {
                    List<volumeType> volType;
                    geom.getVolumeType(testPoints, volType);

                    forAll(insideSurfaces, pointi)
                    {
                        if (insideSurfaces[pointi] == -1)
                        {
                            if (volType[pointi] == volumeType::INSIDE)
                            {
                                insideSurfaces[pointi] = geomi;
                            }
                        }
                    }
                }
            }
        }

        forAll(insideSurfaces, i)
        {
            if (insideSurfaces[i] == -1)
            {
                Pout<< "Flipping cell " << mesh_.cellCentres()[testCells[i]]
                    << endl;
                cuts.flip(testCells[i]);
            }
        }
    }

    // Topo changes container
    polyTopoChange meshMod(mesh_);

    //meshCutAndRemove cutter(mesh_);
    //// Insert mesh refinement into polyTopoChange.
    //cutter.setRefinement
    //(
    //    0,                  //exposedPatci
    //    cuts,
    //    labelList(mesh_.nCells(), 0),    //exposedPatchi
    //    meshMod
    //);


    // Mesh change engine
    Foam::meshCutter cutter(mesh_);
    cutter.setRefinement(cuts, meshMod);

    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);
    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());
    setInstance(mesh_.facesInstance());

    // Update numbering of cells/vertices.
    cutter.updateMesh(map);

    updateMesh(map, identity(mesh_.nFaces()));


    // Detect face on surface:
    // - cut
    // - or all vertices on surface
    if (debug)
    {
        faceSet cutFaces(mesh_, "cutFaces", cutter.addedFaces().size());

        // 1. Get cut faces
        const Map<label>& addedFaces = cutter.addedFaces();
        forAllConstIter(Map<label>, addedFaces, iter)
        {
            cutFaces.insert(iter());
        }

        // 2. Get faces with all vertices snapped
        PackedBoolList isOldVertCut(map().nOldPoints());
        isOldVertCut.set(cutVerts);

        const faceList& faces = mesh_.faces();
        const labelList& pointMap = map().pointMap();

        forAll(faces, facei)
        {
            const face& f = faces[facei];

            bool allSnapped = true;
            forAll(f, fp)
            {
                label oldPointi = pointMap[f[fp]];
                if (!isOldVertCut[oldPointi])
                {
                    allSnapped = false;
                    break;
                }
            }
            if (allSnapped)
            {
                cutFaces.insert(facei);
            }
        }

        Info<< "Writing " << returnReduce(cutFaces.size(), sumOp<label>())
            << " faces to faceSet " << cutFaces.name() << endl;
        cutFaces.write();
    }


    return map;
}


// ************************************************************************* //
