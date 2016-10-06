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
    projectPoints

Description
    Attraction to surface(s)

\*---------------------------------------------------------------------------*/

#include "pointConstraint.H"
#include "argList.H"
#include "fvMesh.H"
#include "searchableSurfaces.H"
#include "plane.H"
#include "OBJstream.H"
#include "lineEdge.H"
#include "lineDivide.H"
#include "gradingDescriptors.H"
#include "triSurfaceMesh.H"
#include "unitConversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Geometric queries
// ~~~~~~~~~~~~~~~~~

void findNearest
(
    const searchableSurfaces& geometry,
    const labelList& surfaces,
    const pointField& start,
    const scalarField& distSqr,
    pointField& near,
    List<pointConstraint>& constraint
)
{
    // Multi-surface findNearest

    vectorField normal;
    List<pointIndexHit> info;

    geometry[surfaces[0]].findNearest(start, distSqr, info);
    geometry[surfaces[0]].getNormal(info, normal);

    // Extract useful info
    near.setSize(info.size());
    forAll(info, i)
    {
        near[i] = info[i].hitPoint();
    }
    constraint.setSize(near.size());

    if (surfaces.size() == 1)
    {
        constraint = pointConstraint();
        forAll(constraint, i)
        {
            constraint[i].applyConstraint(normal[i]);
        }
    }
    else if (surfaces.size() >= 2)
    {
        // Work space
        pointField near1;
        vectorField normal1;

        label surfi = 1;
        for (label iter = 0; iter < 10; iter++)
        {
            constraint = pointConstraint();
            forAll(constraint, i)
            {
                constraint[i].applyConstraint(normal[i]);
            }

            // Find intersection with next surface
            const searchableSurface& s = geometry[surfaces[surfi]];
            s.findNearest(near, distSqr, info);
            s.getNormal(info, normal1);
            near1.setSize(info.size());
            forAll(info, i)
            {
                near1[i] = info[i].hitPoint();
            }

            // Move to intersection
            forAll(near, pointi)
            {
                plane pl0(near[pointi], normal[pointi]);
                plane pl1(near1[pointi], normal1[pointi]);
                plane::ray r(pl0.planeIntersect(pl1));
                vector n = r.dir() / mag(r.dir());

                vector d(r.refPoint()-near[pointi]);
                d -= (d&n)*n;

                near[pointi] += d;
                normal[pointi] = normal1[pointi];

                //Pout<< "point:" << pointi << endl;
                //Pout<< "    pc was:" << constraint[pointi] << endl;
                //Pout<< "    adding:" << normal1[pointi] << endl;
                constraint[pointi].applyConstraint(normal1[pointi]);
                //Pout<< "    pc now:" << constraint[pointi] << endl;
            }

            // Step to next surface
            surfi = surfaces.fcIndex(surfi);
        }
    }
}



// Smoothing
// ~~~~~~~~~

void smooth
(
    const label nIters,
    const edgeList& edges,
    const labelListList& pointEdges,
    const List<pointConstraint>& constraints,
    vectorField& d
)
{
    for (label iter = 0; iter < nIters; iter++)
    {
        pointField oldD(d);

        forAll(pointEdges, pointi)
        {
            const labelList& pEdges = pointEdges[pointi];
            const pointConstraint& pc = constraints[pointi];

            d[pointi] = vector::zero;

            label nNbrs = 0;

            forAll(pEdges, i)
            {
                label nbrPointi = edges[pEdges[i]].otherVertex(pointi);

                if (constraints[nbrPointi].first() >= pc.first())
                {
                    // Feature edges not influenced by normal surface points ...

                    d[pointi] += oldD[nbrPointi];
                    nNbrs++;
                }
            }
            if (nNbrs > 0)
            {
                d[pointi] /= nNbrs;
            }
        }

        d = 0.5*d + 0.5*oldD;
    }
}

void calcFeatureEdges
(
    const triSurface& s,
    const scalar featureAngle,
    PackedBoolList& isBorderEdge
)
{
    isBorderEdge.setSize(s.nEdges());
    isBorderEdge = false;

    scalar cosAngle = Foam::cos(degToRad(featureAngle));

    const labelListList& edgeFaces = s.edgeFaces();
    const vectorField& faceNormals = s.faceNormals();

    forAll(edgeFaces, edgei)
    {
        const labelList& eFaces = edgeFaces[edgei];

        if (eFaces.size() > 2)
        {
            isBorderEdge[edgei] = true;
        }
        else if (eFaces.size() == 2)
        {
            const vector& n0 = faceNormals[eFaces[0]];
            const vector& n1 = faceNormals[eFaces[1]];
            if ((n0&n1) < cosAngle)
            {
                isBorderEdge[edgei] = true;
                Pout<< "MArking feature edge " << edgei
                    << " at " << s.edges()[edgei].centre(s.localPoints())
                    << endl;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    argList::validArgs.append("patch");

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"

    // Read meshing dictionary
    const word dictName("snappyHexMeshDict");
    #include "setSystemMeshDictionaryIO.H"
    const IOdictionary meshDict(dictIO);

    // all surface geometry
    const dictionary& geometryDict = meshDict.subDict("geometry");

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        meshDict.lookupOrDefault("singleRegionName", true)
    );

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const word patchName((IStringStream(args[1])()));
    const polyPatch& pp = pbm[patchName];

    // Find set of patches from the list of regular expressions provided
    //const wordReList patches((IStringStream(args[1])()));
    //const labelList patchIDs(pbm.patchSet(patches).sortedToc());


    if (false)
    {
        OBJstream str("normals.obj");

        pointField patchNear(pp.localPoints());
        List<pointConstraint> patchConstraint(pp.nPoints());
        findNearest
        (
            allGeometry,
            labelList(1, 0),
            pp.localPoints(),
            scalarField(pp.nPoints(), GREAT),
            patchNear,
            patchConstraint
        );
        forAll(patchConstraint, i)
        {
            const vector& n = patchConstraint[i].second();
            str.write(linePointRef(patchNear[i], patchNear[i]+0.1*n));
        }
        return 1;
    }

    if (false)
    {
        // Surfaces to track on intersection
        labelList surfs(2);
        surfs[0] = 0;
        surfs[1] = 1;

        // Create start and end point
        pointField points(2);
        points[0] = point(0, 1, 0);
        points[1] = point(1, 1, 0);

        lineEdge cedge(points, 0, 1);
        gradingDescriptors gds
        (
            gradingDescriptor
            (
                0.5,
                1,
                4
            )
        );

        // Generate a few points
        lineDivide ld(cedge, 20, gds);


        // 1. Find nearest independently
        {
            pointField near;
            List<pointConstraint> constraint;
            findNearest
            (
                allGeometry,
                surfs,
                ld.points(),
                scalarField(ld.points().size(), GREAT),
                near,
                constraint
            );
            OBJstream str("findNearestLine_no_prediction.obj");
            forAll(ld.points(), i)
            {
                str.write(linePointRef(ld.points()[i], near[i]));
            }
        }

        // 2. Track along line
        {
            pointField near(1);
            List<pointConstraint> constraint(1);

            OBJstream str("findNearestLine_with_prediction.obj");

            // Predicted start point
            pointField start(ld.points());
            forAll(start, i)
            {
                Pout<< "inital:" << ld.points()[i] << endl;
                Pout<< "corrected initial:" << start[i] << endl;

                findNearest
                (
                    allGeometry,
                    surfs,
                    pointField(1, start[i]),
                    scalarField(1, GREAT),
                    near,
                    constraint
                );

                str.write(linePointRef(ld.points()[i], start[i]));
                str.write(linePointRef(start[i], near[0]));

                // Do we have a line constraint? If so use the input
                // grading to predict the next starting point
                if (i < start.size()-1 && constraint[0].first() == 2)
                {
                    const vector& edgeN = constraint[0].second();
                    vector d(ld.points()[i+1]-ld.points()[i]);
                    start[i+1] = near[0]+(d&edgeN)*edgeN;
                }
            }
        }

        return 1;
    }



//XXXX
//Problem now is that feature edges are attracted much further than
//internal points and smoothing cannot compensate for that. Need to change order:
//- attraction to feature edges
//- smooth on feature edges and internal
//- calculate new points
//- project internal (with new points)
//- smooth internal
//- calculate new points

    pointField patchNear(pp.localPoints());
    List<pointConstraint> patchConstraint(pp.nPoints());

    scalar relax = 0.1;
    const scalar relaxIncrement = 0.1;
    for (label timeI = 0; timeI < 10; timeI++)
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << endl;

        vectorField displacement;
        {
            patchNear = pp.localPoints();
            patchConstraint = pointConstraint();

            // 1: attract to feature points
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Smooth/interpolate to predict new pp.localPoints()
            // TBD



            // 2: attract to feature edges
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Smooth/interpolate to predict new pp.localPoints()
            {
                labelList surfs(2);
                surfs[0] = 0;
                surfs[1] = 1;
                pointField start(pp.localPoints(), pp.boundaryPoints());

                Pout<< "Attracting " << start.size()
                    << " boundary points to surfaces "
                    << allGeometry[surfs[0]].name() << ' '
                    << allGeometry[surfs[1]].name() << endl;

                pointField boundaryNear;
                List<pointConstraint> boundaryConstraint;
                findNearest
                (
                    allGeometry,
                    surfs,
                    start,
                    scalarField(start.size(), GREAT),
                    boundaryNear,
                    boundaryConstraint
                );

                forAll(pp.boundaryPoints(), i)
                {
                    label pointi = pp.boundaryPoints()[i];
                    patchNear[pointi] = boundaryNear[i];
                    patchConstraint[pointi] = boundaryConstraint[i];
                }
            }

            {
                mkDir(runTime.timePath());
                OBJstream str(runTime.timePath()/"patchNearAfterFeatEdge.obj");
                forAll(patchNear, pointi)
                {
                    const point& pt = patchNear[pointi];
                    str.write(linePointRef(pp.localPoints()[pointi], pt));
                }
            }



            {
                // Smooth
                const vectorField patchDisplacement(patchNear-pp.localPoints());
                displacement = patchDisplacement;
                smooth
                (
                    100,            // nIters,
                    pp.edges(),
                    pp.pointEdges(),
                    patchConstraint,
                    displacement
                );

                {
                    mkDir(runTime.timePath());
                    OBJstream str(runTime.timePath()/"smoothDisp.obj");
                    forAll(displacement, pointi)
                    {
                        const point& pt = pp.localPoints()[pointi];
                        str.write(linePointRef(pt, pt+displacement[pointi]));
                    }
                }

//                // Update displacement to project onto patch
//                forAll(patchConstraint, pointi)
//                {
//                    const pointConstraint& pc = patchConstraint[pointi];
//                    const vector& nearDisp = patchDisplacement[pointi];
//                    vector& d = displacement[pointi];
//
//                    // Displacement d already constrained to surface. Add
//                    // the other direction(s) from the original
//                    // patchDisplacement.
//                    d =
//                        pc.constrainDisplacement(d)
//                      + nearDisp
//                      - pc.constrainDisplacement(nearDisp);
//                }
            }

            {
                mkDir(runTime.timePath());
                OBJstream str(runTime.timePath()/"predictedFeatEdgePoints.obj");
                forAll(displacement, pointi)
                {
                    const point& pt = pp.localPoints()[pointi];
                    str.write(linePointRef(pt, pt+displacement[pointi]));
                }
            }
        }

        // Apply the patch displacement to the mesh points
        {
            pointField newMeshPoints(mesh.points());
            forAll(displacement, pointi)
            {
                newMeshPoints[pp.meshPoints()[pointi]] +=
                    relax*displacement[pointi];
            }
            mesh.movePoints(newMeshPoints);
        }

        runTime.write();

        relax = min(1.0, relax+relaxIncrement);

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }




    relax = 0.1;
    for (label timeI = 0; timeI < 10; timeI++)
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << endl;

        vectorField displacement;
        {
            // 3: attract to single surface
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Determine nearest for all non-boundary points w.r.t. single
            // surface
            {
                const labelList& boundaryPoints = pp.boundaryPoints();

                DynamicList<label> nonBoundaryPoints
                (
                    pp.nPoints()
                   -boundaryPoints.size()
                );
                {
                    PackedBoolList isBoundaryPoint(pp.nPoints());
                    isBoundaryPoint.set(boundaryPoints);

                    forAll(isBoundaryPoint, pointi)
                    {
                        if (!isBoundaryPoint[pointi])
                        {
                            nonBoundaryPoints.append(pointi);
                        }
                    }
                }


                Pout<< "Attracting " << nonBoundaryPoints.size()
                    << " out of " << pp.nPoints()
                    << " patch points to surface " << allGeometry[0].name()
                    << endl;

                pointField nonBoundaryNear;
                List<pointConstraint> nonBoundaryConstraint;
                findNearest
                (
                    allGeometry,
                    labelList(1, 0),
                    pointField(pp.localPoints(), nonBoundaryPoints),
                    scalarField(nonBoundaryPoints.size(), GREAT),
                    nonBoundaryNear,
                    nonBoundaryConstraint
                );
                forAll(nonBoundaryPoints, i)
                {
                    label pointi = nonBoundaryPoints[i];
                    patchNear[pointi] = nonBoundaryNear[i];
                    patchConstraint[pointi] = nonBoundaryConstraint[i];
                }

                {
                    mkDir(runTime.timePath());
                    OBJstream str(runTime.timePath()/"patchNearAfterSurf.obj");
                    forAll(patchNear, pointi)
                    {
                        const point& pt = patchNear[pointi];
                        str.write(linePointRef(pp.localPoints()[pointi], pt));
                    }
                }


                const vectorField patchDisplacement(patchNear-pp.localPoints());
                displacement = patchDisplacement;
                smooth
                (
                    100,            // nIters,
                    pp.edges(),
                    pp.pointEdges(),
                    patchConstraint,
                    displacement
                );

                // Update displacement to project onto patch
                forAll(patchConstraint, pointi)
                {
                    const pointConstraint& pc = patchConstraint[pointi];
                    const vector& nearDisp = patchDisplacement[pointi];
                    vector& d = displacement[pointi];

                    // Displacement d already constrained to surface. Add
                    // the other direction(s) from the original
                    // patchDisplacement.
                    d =
                        pc.constrainDisplacement(d)
                      + nearDisp
                      - pc.constrainDisplacement(nearDisp);
                }
            }
        }


        // Apply the patch displacement to the mesh points
        {
            pointField newMeshPoints(mesh.points());
            forAll(displacement, pointi)
            {
                newMeshPoints[pp.meshPoints()[pointi]] +=
                    relax*displacement[pointi];
            }
            mesh.movePoints(newMeshPoints);
        }

        runTime.write();

        relax = min(1.0, relax+relaxIncrement);

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
