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

tmp<pointField> findNearest(const searchableSurface& s, const pointField& pts)
{
    List<pointIndexHit> info;
    s.findNearest(pts, scalarField(pts.size(), GREAT), info);

    tmp<pointField> tnear(new pointField(pts.size()));
    pointField& near = tnear.ref();
    forAll(info, i)
    {
        near[i] = info[i].hitPoint();
    }
    return tnear;
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

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        vectorField displacement
        (
            findNearest(allGeometry[0], pp.localPoints())
          - pp.localPoints()
        );

        vectorField smoothDisplacement(displacement);
        smooth(pp, smoothDisplacement);

        // TBD: replace tangential component of displacement
        //      with smoothDisplacement
XXXXX


        // Apply the patch displacement to the mesh points
        {
            pointField newMeshPoints(mesh.points());
            forAll(displacement, pointi)
            {
                newMeshPoints[pp.meshPoints()[pointi]] +=
                    0.5*displacement[pointi];
            }
            mesh.movePoints(newMeshPoints);
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
