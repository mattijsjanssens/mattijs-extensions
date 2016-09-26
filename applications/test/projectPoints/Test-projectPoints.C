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
    writePatchDistance

Description
    Write topological distance to wall

\*---------------------------------------------------------------------------*/

#include "pointConstraint.H"
#include "argList.H"
#include "fvMesh.H"
#include "searchableSurfaces.H"
#include "plane.H"
#include "OBJstream.H"

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

                Pout<< "point:" << pointi << endl;
                Pout<< "    pc was:" << constraint[pointi] << endl;
                Pout<< "    adding:" << normal1[pointi] << endl;
                constraint[pointi].applyConstraint(normal1[pointi]);
                Pout<< "    pc now:" << constraint[pointi] << endl;
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
                d[pointi] /= pEdges.size();
                d[pointi] = pc.constrainDisplacement(d[pointi]);
            }
        }

        d = 0.5*d + 0.5*oldD;
    }
}


//XXXX
//void findNearest
//(
//    const searchableSurface& s,
//    const pointField& pts,
//    pointField& near,
//    vectorField& normal
//)
//{
//    List<pointIndexHit> info;
//    s.findNearest(pts, scalarField(pts.size(), GREAT), info);
//    s.getNormal(info, normal);
//
//    near.setSize(info.size());
//    forAll(info, i)
//    {
//        near[i] = info[i].hitPoint();
//    }
//}
//void smooth(const primitivePatch& p, vectorField& d)
//{
//    const labelListList& pointEdges = p.pointEdges();
//    const edgeList& edges = p.edges();
//
//    for (label iter = 0; iter < 100; iter++)
//    {
//        pointField oldD(d);
//
//        forAll(pointEdges, pointi)
//        {
//            const labelList& pEdges = pointEdges[pointi];
//
//            d[pointi] = vector::zero;
//            forAll(pEdges, i)
//            {
//                label nbrPointi = edges[pEdges[i]].otherVertex(pointi);
//
//                d[pointi] += oldD[nbrPointi];
//            }
//            d[pointi] /= pEdges.size();
//        }
//
//        d = 0.5*d + 0.5*oldD;
//    }
//}
//void smoothAttraction
//(
//    const searchableSurface& s,
//    const primitivePatch& pp,
//    vectorField& displacement
//)
//{
//    pointField surfNearest;
//    vectorField surfNormal;
//    findNearest(s, pp.localPoints(), surfNearest, surfNormal);
//
//    vectorField surfDisplacement(surfNearest-pp.localPoints());
//
//    displacement = surfDisplacement;
//    smooth(pp, displacement);
//
//    // TBD: replace normal component of displacement with that
//    //      of surfDisplacement
//    displacement -= (displacement & surfNormal) * surfNormal;
//    displacement += (surfDisplacement & surfNormal) * surfNormal;
//}
//XXXXX

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


//    {
//        labelList surfs(2);
//        surfs[0] = 0;
//        surfs[1] = 1;
//
//        // Detect points on outside of patch and attract them to the
//        // nearest feature lines
//        pointField start(pp.localPoints(), pp.boundaryPoints());
//
//        pointField near;
//        List<pointConstraint> constraint;
//        findNearest
//        (
//            allGeometry,
//            surfs,
//            start,
//            scalarField(start.size(), GREAT),
//            near,
//            constraint
//        );
//
//        DebugVar(constraint);
//
//        OBJstream str("findNearestLine.obj");
//        forAll(start, i)
//        {
//            str.write(linePointRef(start[i], near[i]));
//        }
//    }



//    // Scan patches to find feature line points and feature points:
//    // construct map from mesh point to set of patches
//    Map<labelList> pointToPatches;
//    forAll(patchIDs, i)
//    {
//        label patchi = patchIDs[i];
//        const labelList& meshPPoints = pbm[patchi].meshPoints();
//        const labelList& boundaryPoints = pbm[patchi].boundaryPoints();
//        forAll(boundaryPoints, bPointi)
//        {
//            label meshPointi = meshPPoints[bPointi];
//            Map<labelList>::iterator iter = pointToPatches.find(meshPointi);
//            if (iter == pointToPatches.end())
//            {
//                pointToPatches.insert(meshPointi, labelList(1, i));
//            }
//            else
//            {
//                iter().append(i);
//            }
//        }
//    }

    scalar relax = 0.1;
    scalar relaxIncrement = 0.1;


XXXX
Problem now is that feature edges are attracted much further than
internal points and smoothing cannot compensate for that. Need to change order:
- attraction to feature edges
- smooth on feature edges and internal
- calculate new points
- project internal (with new points)
- smooth internal
- calculate new points


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        vectorField displacement;
        {
//          smoothAttraction(allGeometry[0], pp, displacement);

            // Determine nearest for all points w.r.t. single surface
            Pout<< "Attracting " << pp.nPoints()
                << " boundary points to surfae " << allGeometry[0].name()
                << endl;
            pointField patchNear;
            List<pointConstraint> patchConstraint;
            findNearest
            (
                allGeometry,
                labelList(1, 0),
                pp.localPoints(),
                scalarField(pp.nPoints(), GREAT),
                patchNear,
                patchConstraint
            );

            // Determine nearest for boundary points w.r.t. two surfaces and
            // override
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

            forAll(patchConstraint, pointi)
            {
                const pointConstraint& pc = patchConstraint[pointi];
                const vector& nearDisp = patchDisplacement[pointi];
                vector& d = displacement[pointi];

                // Displacement d already constrained to surface. Add
                // the other direction(s) from the original patchDisplacement.

                d += nearDisp - pc.constrainDisplacement(nearDisp);
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
