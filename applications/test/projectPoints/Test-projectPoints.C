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


// Two surface findNearest
void findNearest
(
    const searchableSurface& s0,
    const searchableSurface& s1,
    const point& start,
    point& near
//,
//    vector& normal0,
//    vector& normal1
)
{
    point pt(start);
    vector n;
    {
        List<pointIndexHit> info0;
        s0.findNearest(pointField(1, pt), scalarField(1, GREAT), info0);
        vectorField normal1;
        s.getNormal(info0, normal0);

        pt = info0[0].hitPoint();
        n = normal0[0];
    }

    // Find intersection with other surface
    point pt1;
    vector n1;
    {
        List<pointIndexHit> info1;
        s1.findNearest(pointField(1, pt), scalarField(1, GREAT), info1);
        vectorField normal1;
        s.getNormal(info1, normal1);

        // Move to intersection
        plane pl0(pt, n);
        plane pl1(info1[0].hitPoint(), normal1[0]);
        plane::ray r(pl0.planeIntersect(pl1));
        vector n = r.dir() / mag(r.dir());

        vector d(r.refPoint()-pt);
        d -= (d&n)*n;
        pt1 = pt + d;

        n1 = normal1[0];
    }


    near.setSize(info.size());
    forAll(info, i)
    {
        near[i] = info[i].hitPoint();
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
