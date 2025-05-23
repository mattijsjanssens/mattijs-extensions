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

Application
    cutMeshWithSurface

Description
    cut mesh at the intersections with a surface. Either cut cells
    into two (uncomment meshCut) or cut-and-remove the neighbour side
    of the cut (uncomment meshCutAndRemove)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "faceSet.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "triSurfaceMesh.H"
#include "DynamicField.H"

#include "polyTopoChange.H"
#include "cellCuts.H"
//#include "meshCutAndRemove.H"
#include "meshCutter.H"
#include "pointSet.H"
#include "cellSet.H"
#include "OBJstream.H"
#include "zeroGradientFvPatchFields.H"
#include "fvMeshSubset.H"
#include "regionSplit.H"
#include "surfaceFeatures.H"
#include "treeDataFace.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

PackedBoolList isFaceCut(const polyMesh& mesh, const triSurfaceMesh& surfMesh)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const pointField& cellCentres = mesh.cellCentres();

    pointField start(mesh.nFaces());
    pointField end(mesh.nFaces());
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        start[facei] = cellCentres[own[facei]];
        end[facei] = cellCentres[nei[facei]];
    }
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        start[facei] = cellCentres[own[facei]];
        vector d(mesh.faceCentres()[facei]-start[facei]);
        end[facei] = start[facei]+2*d;
    }

    List<pointIndexHit> info;
    surfMesh.findLineAny(start, end, info);

    PackedBoolList isCut(mesh.nFaces());
    forAll(info, facei)
    {
        if (info[facei].hit())
        {
            isCut[facei] = true;
        }
    }
    return isCut;
}


label createBaffle
(
    const fvMesh& mesh,
    const label facei,
    const label ownPatch,
    const label neiPatch,
    polyTopoChange& meshMod
)
{
    const face& f = mesh.faces()[facei];
    label zoneID = mesh.faceZones().whichZone(facei);
    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];
        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }

    meshMod.modifyFace
    (
        f,                          // modified face
        facei,                      // label of face
        mesh.faceOwner()[facei],    // owner
        -1,                         // neighbour
        false,                      // face flip
        ownPatch,                   // patch for face
        zoneID,                     // zone for face
        zoneFlip                    // face flip in zone
    );


    label dupFacei = -1;

    if (mesh.isInternalFace(facei))
    {
        if (neiPatch == -1)
        {
            FatalErrorInFunction
                << "No neighbour patch for internal face " << facei
                << " fc:" << mesh.faceCentres()[facei]
                << " ownPatch:" << ownPatch << abort(FatalError);
        }

        bool reverseFlip = false;
        if (zoneID >= 0)
        {
            reverseFlip = !zoneFlip;
        }

        dupFacei = meshMod.addFace
        (
            f.reverseFace(),            // modified face
            mesh.faceNeighbour()[facei],// owner
            -1,                         // neighbour
            -1,                         // masterPointID
            -1,                         // masterEdgeID
            facei,                      // masterFaceID,
            true,                       // face flip
            neiPatch,                   // patch for face
            zoneID,                     // zone for face
            reverseFlip                 // face flip in zone
        );
    }
    return dupFacei;
}


// Do feature points
void snapToFeaturePoints
(
    fvMesh& mesh,
    const triSurfaceMesh& surfMesh,
    const surfaceFeatures& feats,
    PackedBoolList& isCutVert
)
{
    treeBoundBox bb(surfMesh.bounds());
    bb.inflate(1e-3);

    // Find cells containing the feature points
    const labelList& featurePoints = feats.featurePoints();

    Info<< "For surface " << surfMesh.searchableSurface::name()
        << " detected " << featurePoints.size() << " feature points." << endl;

    pointField newPoints(mesh.points());

    forAll(featurePoints, feati)
    {
        const point& featPt = surfMesh.localPoints()[featurePoints[feati]];

        label celli = mesh.findCell(featPt);
        if (celli != -1)
        {
            const labelList& cPoints = mesh.cellPoints(celli);


            scalar minDistSqr = GREAT;
            label minFp = -1;
            forAll(cPoints, fp)
            {
                label pointi = cPoints[fp];
                const point& pt = mesh.points()[pointi];
                scalar distSqr = magSqr(pt-featPt);
                if (distSqr < minDistSqr)
                {
                    minDistSqr = distSqr;
                    minFp = fp;
                }
            }
            if (minFp != -1)
            {
                label pointi = cPoints[minFp];
                newPoints[pointi] = featPt;
                isCutVert[pointi] = true;
            }
        }
    }
    mesh.movePoints(newPoints);
    const_cast<Time&>(mesh.time())++;
    Info<< "Writing feature-point snapped mesh to time "
        << mesh.time().timeName() << endl;
    mesh.write();
}


// Do feature edges
bool onCell(const fvMesh& mesh, const label facei, const label otherFacei)
{
    label own = mesh.faceOwner()[facei];
    label otherOwn = mesh.faceOwner()[otherFacei];

    if (otherOwn == own)
    {
        return true;
    }

    if (mesh.isInternalFace(facei))
    {
        label nei = mesh.faceNeighbour()[facei];
        if (otherOwn == nei)
        {
            return true;
        }

        if (mesh.isInternalFace(otherFacei))
        {
            label otherNei = mesh.faceNeighbour()[otherFacei];
            if (otherNei == own || otherNei == nei)
            {
                return true;
            }
        }
        return false;
    }
    else
    {
        if (mesh.isInternalFace(otherFacei))
        {
            label otherNei = mesh.faceNeighbour()[otherFacei];
            if (otherNei == own)
            {
                return true;
            }
        }
        return false;
    }
}
void snapToFeatureEdges
(
    fvMesh& mesh,
    const triSurfaceMesh& surfMesh,
    const surfaceFeatures& feats,
    PackedBoolList& isCutVert,
    labelList& whichFace
)
{
    treeBoundBox bb(surfMesh.bounds());
    bb.inflate(1e-3);

    indexedOctree<treeDataFace> tree
    (
        treeDataFace
        (
            false,                  // do not cache bb
            mesh,
            identity(mesh.nFaces())
        ),
        bb,     // overall search domain
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );

    // Find intersections of feature edges with mesh faces
    const labelList& featureEdges = feats.featureEdges();

    Info<< "For surface " << surfMesh.searchableSurface::name()
        << " detected " << featureEdges.size() << " feature edges." << endl;

    {
        OBJstream str(mesh.time().path()/"featureEdges.obj");
        forAll(featureEdges, feati)
        {
            const edge& surfE = surfMesh.edges()[featureEdges[feati]];
            const point& start = surfMesh.localPoints()[surfE[0]];
            const point& end = surfMesh.localPoints()[surfE[1]];
            str.write(linePointRef(start, end));
        }
    }



    pointField newPoints(mesh.points());
    scalarField snapDistSqr(mesh.nPoints(), GREAT);

    // Index of face that caused feature attraction
    whichFace.setSize(mesh.nPoints());
    whichFace = -1;

    OBJstream str(mesh.time().path()/"featureIntersections.obj");
    forAll(featureEdges, feati)
    {
        label surfEdgei = featureEdges[feati];
        const edge& surfE = surfMesh.edges()[surfEdgei];
        point start = surfMesh.localPoints()[surfE[0]];
        const point& end = surfMesh.localPoints()[surfE[1]];
        const vector vec(end-start);
        const vector smallVec(1e-6*(vec));

        while (magSqr(end-start) > SMALL)
        {
            pointIndexHit info = tree.findLine(start, end-smallVec);

            if (!info.hit())
            {
                break;
            }

            str.write(info.hitPoint());

            // Find the nearest vertex to attract
            const label facei = info.index();
            const face& f = mesh.faces()[facei];

            scalar minDistSqr = GREAT;
            label minFp = -1;

            forAll(f, fp)
            {
                label pointi = f[fp];
                const point& pt = mesh.points()[pointi];
                scalar distSqr = magSqr(pt-info.hitPoint());
                if (distSqr < minDistSqr)
                {
                    minDistSqr = distSqr;
                    minFp = fp;
                }
            }
            if (minFp != -1)
            {
                label pointi = f[minFp];
                if (!isCutVert[pointi] && (minDistSqr < snapDistSqr[pointi]))
                {
                    snapDistSqr[pointi] = minDistSqr;
                    newPoints[pointi] = info.hitPoint();

                    isCutVert[pointi] = true;

                    if (whichFace[pointi] == -1)
                    {
                        whichFace[pointi] = facei;
                    }
                    else
                    {
                        whichFace[pointi] = -2;
                    }
                }
            }

            start = info.hitPoint() + smallVec;
        }
    }


    // Shift whole planes
    {
        OBJstream str(mesh.time().path()/"planeSnap.obj");

        forAll(isCutVert, pointi)
        {
            if (isCutVert[pointi] && whichFace[pointi] >= 0)
            {
                label facei = whichFace[pointi];
                //label own = mesh.faceOwner()[facei];
                //label nei = mesh.faceNeighbour()[facei];

                const labelList& pFaces = mesh.pointFaces()[pointi];
                forAll(pFaces, pFacei)
                {
                    label otherFacei = pFaces[pFacei];

                    if (otherFacei != facei && onCell(mesh, facei, otherFacei))
                    {
                        label nAttract = 0;
                        const face& otherF = mesh.faces()[otherFacei];
                        forAll(otherF, fp)
                        {
                            if (isCutVert[otherF[fp]])
                            {
                                nAttract++;
                            }
                        }

                        if (nAttract == 1)
                        {
                            Pout<< "found face:"
                                << mesh.faceCentres()[otherFacei];
                            str.write(otherF, mesh.points(), false);

                            const vector d
                            (
                                newPoints[pointi]
                               -mesh.points()[pointi]
                            );

                            forAll(otherF, fp)
                            {
                                if (otherF[fp] != pointi)
                                {
                                    newPoints[otherF[fp]] += d;
                                }
                            }
                        }
                    }
                }
            }
        }
    }



    mesh.movePoints(newPoints);
    const_cast<Time&>(mesh.time())++;
    Info<< "Writing feature-edge snapped mesh to time "
        << mesh.time().timeName() << endl;
    mesh.write();
}


//void planeSnap
//(
//    fvMesh& mesh,
//    const triSurfaceMesh& surfMesh,
//    const surfaceFeatures& feats,
//    PackedBoolList& isCutVert,
//    labelList& whichFace
//)
//{
//    OBJstream str(mesh.time().path()/"planeSnap.obj");
//
//    forAll(isCutVert, pointi)
//    {
//        if (isCutVert[pointi] && whichFace[pointi] >= 0)
//        {
//            label facei = whichFace[pointi];
//            //label own = mesh.faceOwner()[facei];
//            //label nei = mesh.faceNeighbour()[facei];
//
//            const labelList& pFaces = mesh.pointFaces()[pointi];
//            forAll(pFaces, pFacei)
//            {
//                label otherFacei = pFaces[pFacei];
//
//                if (otherFacei != facei && onCell(mesh, facei, otherFacei))
//                {
//                    label nAttract = 0;
//                    const face& otherF = mesh.faces()[otherFacei];
//                    forAll(otherF, fp)
//                    {
//                        if (isCutVert[otherF[fp]])
//                        {
//                            nAttract++;
//                        }
//                    }
//
//                    if (nAttract == 1)
//                    {
//                        Pout<< "found face:" << mesh.faceCentres()[otherFacei];
//                        str.write(otherF, mesh.points(), false);
//                    }
//                }
//            }
//        }
//    }
//}


int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    argList::validArgs.append("triSurfaceMesh");

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const fileName surfFileName = args[1];
    const bool overwrite = args.optionFound("overwrite");



    triSurfaceMesh surfMesh
    (
        IOobject
        (
            surfFileName,         // name
            runTime.constant(),   // instance
            "triSurface",         // local
            runTime,              // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
        //45.0
    );


    // Find the edges intersecting the surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<label> cutVerts;
    DynamicList<label> cutEdges;
    DynamicField<scalar> cutEdgeWeights;


    // Option 1: always cut edges
    if (false)
    {
        const edgeList& edges = mesh.edges();
        const pointField& points = mesh.points();

        List<pointIndexHit> info;
        {
            vectorField start(edges.size());
            vectorField end(edges.size());
            forAll(edges, edgei)
            {
                const edge& e = edges[edgei];
                start[edgei] = points[e[0]];
                end[edgei] = points[e[1]];
            }
            surfMesh.findLineAny(start, end, info);
            //vectorField normal;
            //surfMesh.getNormal(info, normal);
        }


        forAll(info, edgei)
        {
            if (info[edgei].hit())
            {
                const edge& e = edges[edgei];
                const point& pt = info[edgei].hitPoint();

                vector eVec(e.vec(points));
                scalar f = eVec&(pt-points[e.start()])/magSqr(eVec);

                // if (f < 0.1)
                // {
                //     cutVerts.append(e.start());
                // }
                // else if (f > 0.9)
                // {
                //     cutVerts.append(e.end());
                // }
                // else
                {
                    cutEdges.append(edgei);
                    cutEdgeWeights.append(f);
                }
            }
        }
    }

    // Option 2: approximate the surface through motion and cutting
    if (false)
    {
        const labelList& own = mesh.faceOwner();
        const labelList& nei = mesh.faceNeighbour();
        const pointField& cellCentres = mesh.cellCentres();

        pointField start(mesh.nFaces());
        pointField end(mesh.nFaces());
        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            start[facei] = cellCentres[own[facei]];
            end[facei] = cellCentres[nei[facei]];
        }
        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            start[facei] = cellCentres[own[facei]];
            vector d(mesh.faceCentres()[facei]-start[facei]);
            end[facei] = start[facei]+2*d;
        }

        List<pointIndexHit> info;
        surfMesh.findLineAny(start, end, info);
        vectorField normal;
        surfMesh.getNormal(info, normal);
    }


    // Option 3: always cut vertices by moving edges along
    {
        const edgeList& edges = mesh.edges();
        const pointField& points = mesh.points();

        PackedBoolList isCutVert(points.size());


        surfaceFeatures feats(surfMesh, 180.0-45.0);

        snapToFeaturePoints(mesh, surfMesh, feats, isCutVert);
        labelList whichFace;
        snapToFeatureEdges(mesh, surfMesh, feats, isCutVert, whichFace);

        //{
        //    planeSnap
        //    (
        //        mesh,
        //        surfMesh,
        //        feats,
        //        isCutVert,
        //        whichFace
        //    );
        //}


        List<pointIndexHit> info;
        vectorField normal;
        for (label iter = 0;; iter++)
        {
            // Update intersections
            {
                vectorField start(edges.size());
                vectorField end(edges.size());
                forAll(edges, edgei)
                {
                    const edge& e = edges[edgei];
                    start[edgei] = points[e[0]];
                    end[edgei] = points[e[1]];
                }
                surfMesh.findLineAny(start, end, info);
                surfMesh.getNormal(info, normal);
            }

            if (iter == 10)
            {
                break;
            }


            pointField newPoints(points);
            label nAdjusted = 0;

            const labelListList& pointEdges = mesh.pointEdges();
            forAll(pointEdges, pointi)
            {
                if (!isCutVert[pointi])
                {
                    const point& pt = points[pointi];
                    const labelList& pEdges = pointEdges[pointi];

                    // Get the nearest intersection
                    label minEdgei = -1;
                    scalar minFraction = 0.5;   // Harpoon 0.25; // Samm?
                    forAll(pEdges, pEdgei)
                    {
                        label edgei = pEdges[pEdgei];
                        if (info[edgei].hit())
                        {
                            const point& hitPt = info[edgei].hitPoint();

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
                        newPoints[pointi] = info[minEdgei].hitPoint();
                        nAdjusted++;
                    }
                }
            }

            Info<< "Iter:" << iter << " adjusting " << nAdjusted << " points"
                << endl;

            if (nAdjusted == 0)
            {
                break;
            }
            mesh.movePoints(newPoints);


            if (!overwrite)
            {
                runTime++;
            }
            else
            {
                mesh.setInstance(oldInstance);
            }

            // Take over refinement levels and write to new time directory.
            Info<< "Writing mesh to time " << runTime.timeName() << endl;
            mesh.write();
        }


        // Filter out edge cuts if endpoints have been cut
        cutVerts.clear();
        cutEdges.clear();
        cutEdgeWeights.clear();
        forAll(info, edgei)
        {
            if (info[edgei].hit())
            {
                const edge& e = edges[edgei];

                if (!isCutVert[e[0]] && !isCutVert[e[1]])
                {
                    const point& pt = info[edgei].hitPoint();

                    vector eVec(e.vec(points));
                    scalar f = eVec&(pt-points[e[0]])/magSqr(eVec);

                    if (f < 0.01)
                    {
                        isCutVert[e[0]] = true;
                    }
                    else if (f > 0.99)
                    {
                        isCutVert[e[1]] = true;
                    }
                }
            }
        }

        forAll(info, edgei)
        {
            if (info[edgei].hit())
            {
                const edge& e = edges[edgei];

                if (!isCutVert[e[0]] && !isCutVert[e[1]])
                {
                    const point& pt = info[edgei].hitPoint();
                    vector eVec(e.vec(points));
                    scalar f = eVec&(pt-points[e[0]])/magSqr(eVec);

                    cutEdges.append(edgei);
                    cutEdgeWeights.append(f);
                }
            }
        }

        forAll(isCutVert, pointi)
        {
            if (isCutVert[pointi])
            {
                cutVerts.append(pointi);
            }
        }
    }


    {
        const edgeList& edges = mesh.edges();
        const pointField& points = mesh.points();

        OBJstream str(runTime.path()/"allCuts.obj");
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
    cellCuts cuts(mesh, cutVerts, cutEdges, cutEdgeWeights);

    DebugVar(cuts.nLoops());

    // See if we need to flip anything
    {
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
                testPoints.append(mesh.points()[anchors[0]]);
            }
        }

        List<volumeType> side;
        surfMesh.getVolumeType(testPoints, side);
        forAll(side, i)
        {
            if (side[i] == volumeType::OUTSIDE)
            {
                Pout<< "Flipping cell " << mesh.cellCentres()[testCells[i]]
                    << endl;
                cuts.flip(testCells[i]);
            }
        }
    }

    // Dump failed cuts to cellSet
    {
        cellSet candidates(mesh, "failedCandidates", mesh.nCells()/100);
        forAll(cutEdges, i)
        {
            const labelList& cEdges = mesh.edgeCells()[cutEdges[i]];
            candidates.insert(cEdges);
        }
        DebugVar(candidates.size());
        forAll(cutVerts, i)
        {
            const labelList& pCells = mesh.pointCells()[cutVerts[i]];
            candidates.insert(pCells);
        }
        DebugVar(candidates.size());
        const labelListList& cellLoops = cuts.cellLoops();
        forAll(cellLoops, celli)
        {
            if (cellLoops[celli].size())
            {
                candidates.erase(celli);
            }
        }
        DebugVar(candidates.size());

        candidates.write();
    }


//     // Read objects in time directory
//     IOobjectList objects(mesh, runTime.timeName());
//
//     // Read vol fields.
//     PtrList<volScalarField> vsFlds;
//     ReadFields(mesh, objects, vsFlds);
//
//     PtrList<volVectorField> vvFlds;
//     ReadFields(mesh, objects, vvFlds);
//
//     PtrList<volSphericalTensorField> vstFlds;
//     ReadFields(mesh, objects, vstFlds);
//
//     PtrList<volSymmTensorField> vsymtFlds;
//     ReadFields(mesh, objects, vsymtFlds);
//
//     PtrList<volTensorField> vtFlds;
//     ReadFields(mesh, objects, vtFlds);
//
//- Reads meshPhi which does not get mapped sufficiently (patch does
//  not get extended for added faces)
//    // Read surface fields.
//    PtrList<surfaceScalarField> ssFlds;
//    ReadFields(mesh, objects, ssFlds);
//
//    PtrList<surfaceVectorField> svFlds;
//    ReadFields(mesh, objects, svFlds);
//
//    PtrList<surfaceSphericalTensorField> sstFlds;
//    ReadFields(mesh, objects, sstFlds);
//
//    PtrList<surfaceSymmTensorField> ssymtFlds;
//    ReadFields(mesh, objects, ssymtFlds);
//
//    PtrList<surfaceTensorField> stFlds;
//    ReadFields(mesh, objects, stFlds);



    // Topo changes container
    polyTopoChange meshMod(mesh);

    //meshCutAndRemove cutter(mesh);
    //// Insert mesh refinement into polyTopoChange.
    //cutter.setRefinement
    //(
    //    0,                  //exposedPatci
    //    cuts,
    //    labelList(mesh.nCells(), 0),    //exposedPatchi
    //    meshMod
    //);


    // Mesh change engine
    meshCutter cutter(mesh);
    cutter.setRefinement(cuts, meshMod);



    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    mesh.updateMesh(morphMap);

    // Move mesh (since morphing does not do this)
    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    // Update numbering of cells/vertices.
    cutter.updateMesh(morphMap);

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    // Take over refinement levels and write to new time directory.
    Info<< "Writing mesh to time " << runTime.timeName() << endl;
    mesh.write();


    // Detect face on surface:
    // - cut
    // - or all vertices on surface
    faceSet cutFaces(mesh, "cutFaces", cutter.addedFaces().size());
    {
        // 1. Get cut faces
        const Map<label>& addedFaces = cutter.addedFaces();
        forAllConstIter(Map<label>, addedFaces, iter)
        {
            cutFaces.insert(iter());
        }

        // 2. Get faces with all vertices snapped
        PackedBoolList isOldVertCut(morphMap().nOldPoints());
        isOldVertCut.set(cutVerts);

        const faceList& faces = mesh.faces();
        const labelList& pointMap = morphMap().pointMap();

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

    {
        volScalarField fld
        (
            IOobject
            (
                "cellMap",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(fld, celli)
        {
            fld[celli] = 1.0*morphMap().cellMap()[celli];
        }
        fld.write();
    }


    // Create some baffles
    {
        // Topo changes container
        polyTopoChange meshMod(mesh);

        forAllConstIter(faceSet, cutFaces, iter)
        {
            createBaffle
            (
                mesh,
                iter.key(),
                0,          // ownPatch,
                0,          // neiPatch,
                meshMod
            );
        }

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        mesh.updateMesh(morphMap);

        // Move mesh (since morphing does not do this)
        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        // Update numbering of cells/vertices.
        cutter.updateMesh(morphMap);

        if (!overwrite)
        {
            runTime++;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        // Take over refinement levels and write to new time directory.
        Info<< "Writing baffle mesh to time " << runTime.timeName() << endl;
        mesh.write();
    }


    // Split into separate meshes
    {
        regionSplit cellRegion(mesh);
        Info<< "Mesh has " << cellRegion.nRegions() << " regions" << endl;

        fvMeshSubset subsetter(mesh);

        for (label regioni = 0; regioni < cellRegion.nRegions(); regioni++)
        {
            subsetter.setLargeCellSubset(cellRegion, regioni, -1);

            runTime++;

            Info<< "Writing region " << regioni
                << " mesh to time " << runTime.timeName() << endl;
            subsetter.subMesh().write();
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
