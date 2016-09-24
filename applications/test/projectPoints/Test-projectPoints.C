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

#include "argList.H"
#include "fvMesh.H"
#include "searchableSurfaces.H"
#include "plane.H"
#include "OBJstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//tmp<pointField> findNearest(const searchableSurface& s, const pointField& pts)
//{
//    List<pointIndexHit> info;
//    s.findNearest(pts, scalarField(pts.size(), GREAT), info);
//
//    tmp<pointField> tnear(new pointField(pts.size()));
//    pointField& near = tnear.ref();
//    forAll(info, i)
//    {
//        near[i] = info[i].hitPoint();
//    }
//    return tnear;
//}
void findNearest
(
    const searchableSurface& s,
    const pointField& pts,
    pointField& near,
    vectorField& normal
)
{
    List<pointIndexHit> info;
    s.findNearest(pts, scalarField(pts.size(), GREAT), info);
    s.getNormal(info, normal);

    near.setSize(info.size());
    forAll(info, i)
    {
        near[i] = info[i].hitPoint();
    }
}
void smooth(const primitivePatch& p, vectorField& d)
{
    const labelListList& pointEdges = p.pointEdges();
    const edgeList& edges = p.edges();

    for (label iter = 0; iter < 100; iter++)
    {
        pointField oldD(d);

        forAll(pointEdges, pointi)
        {
            const labelList& pEdges = pointEdges[pointi];

            d[pointi] = vector::zero;
            forAll(pEdges, i)
            {
                label nbrPointi = edges[pEdges[i]].otherVertex(pointi);

                d[pointi] += oldD[nbrPointi];
            }
            d[pointi] /= pEdges.size();
        }

        d = 0.5*d + 0.5*oldD;
    }
}
void smoothAttraction
(
    const searchableSurface& s,
    const primitivePatch& pp,
    vectorField& displacement
)
{
    pointField surfNearest;
    vectorField surfNormal;
    findNearest(s, pp.localPoints(), surfNearest, surfNormal);

    vectorField surfDisplacement(surfNearest-pp.localPoints());

    displacement = surfDisplacement;
    smooth(pp, displacement);

    // TBD: replace normal component of displacement with that
    //      of surfDisplacement
    displacement -= (displacement & surfNormal) * surfNormal;
    displacement += (surfDisplacement & surfNormal) * surfNormal;
}


// one or two surface findNearest
void findNearest
(
    const searchableSurfaces& geometry,
    const labelList& surfaces,
    const pointField& start,
    pointField& near,
    vectorField& normal
)
{
    findNearest(geometry[surfaces[0]], start, near, normal);

    // Work space
    pointField near1;
    vectorField normal1;

    if (surfaces.size() == 2)
    {
        label surfi = 1;
        for (label iter = 0; iter < 10; iter++)
        {
            // Find intersection with next surface
            findNearest(geometry[surfaces[surfi]], near, near1, normal1);

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
            }

            // Step to next surface
            surfi = surfaces.fcIndex(surfi);
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

    const word patchName((IStringStream(args[1])()));

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const polyPatch& pp = pbm[patchName];


    {
        labelList surfs(2);
        surfs[0] = 0;
        surfs[1] = 1;

        // Detect points on outside of patch and attract them to the
        // nearest feature lines
        pointField start(pp.localPoints(), pp.boundaryPoints());

        pointField near;
        vectorField normal;
        findNearest
        (
            allGeometry,
            surfs,
            start,
            near,
            normal
        );

        OBJstream str("findNearest.obj");
        forAll(start, i)
        {
            str.write(linePointRef(start[i], near[i]));
        }

        return 1;
    }





    scalar relax = 0.1;
    scalar relaxIncrement = 0.1;


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        vectorField displacement;
        smoothAttraction(allGeometry[0], pp, displacement);

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
