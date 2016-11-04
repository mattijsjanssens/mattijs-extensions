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

Description
    writeSurfaceDistance

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "pointMesh.H"
#include "OSspecific.H"
#include "IFstream.H"
//#include "pointEdgePoint.H"
#include "pointData.H"
#include "PointEdgeWave.H"
#include "smoothTriSurfaceMesh.H"
#include "OBJstream.H"
#include "DynamicField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("triSurfaceMesh");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const fileName surfFileName = args[1];

    smoothTriSurfaceMesh surfMesh
    (
        IOobject
        (
            surfFileName,         // name
            runTime.constant(),   // instance
            "triSurface",         // local
            runTime,              // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        45.0
    );


    // Find the edge intersecting the surface
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    vectorField start(edges.size());
    vectorField end(edges.size());
    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        start[edgei] = points[e[0]];
        end[edgei] = points[e[1]];
    }
    List<pointIndexHit> info;
    surfMesh.findLineAny(start, end, info);
    vectorField normal;
    surfMesh.getNormal(info, normal);

    
    // Seed distance from intersected edges
    DynamicList<pointData> wallInfo;
    DynamicList<label> wallPoints;

    forAll(info, edgei)
    {
        if (info[edgei].hit())
        {
            const point& pt = info[edgei].hitPoint();
            //const label triI = info[edgei].index();
            //const point& fc = surfMesh.faceCentres()[triI];

            const edge& e = edges[edgei];

            wallPoints.append(e[0]);
            wallInfo.append
            (
                pointData
                (
                    pt, //start[edgei],
                    magSqr(start[edgei]-pt),
                    scalar(info[edgei].index()),    //scalar(e[0]),
                    normal[edgei]
                )
            );

            wallPoints.append(e[1]);
            wallInfo.append
            (
                pointData
                (
                    pt, //end[edgei],
                    magSqr(end[edgei]-pt),
                    scalar(info[edgei].index()),    //scalar(e[1]),
                    normal[edgei]
                )
            );
        }
    }

    Info<< "Seeded " << returnReduce(wallPoints.size(), sumOp<label>())
        << " points" << nl << endl;



    // Current info on points
    List<pointData> allPointInfo(mesh.nPoints());

    // Current info on edges
    List<pointData> allEdgeInfo(mesh.nEdges());

    PointEdgeWave<pointData> wallCalc
    (
        mesh,
        wallPoints,
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh.nPoints()  // max iterations
    );


    pointScalarField psf
    (
        IOobject
        (
            "wallDist",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedScalar("wallDist", dimLength, 0.0)
    );

    forAll(allPointInfo, pointI)
    {
        psf[pointI] = Foam::sqrt(allPointInfo[pointI].distSqr());
    }

    Info<< "Writing wallDist pointScalarField to " << runTime.value()
        << endl;

    psf.write();


    // Extract new surface points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicField<point> isoPoints;
    DynamicField<point> surfPoints;
    DynamicField<vector> surfNormals;



    // Plot lines from iso at 50 back to originating surface
    OBJstream str(runTime.path()/"surface_to_iso.obj");


    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        const point& origin0 = allPointInfo[e[0]].origin();
        const point& origin1 = allPointInfo[e[1]].origin();

        scalar delta = psf[e[1]]-psf[e[0]];

        if (mag(delta) > SMALL)
        {
            // Note: small perturb since some points at exactly 50.
            scalar f = (49.937-psf[e[0]])/delta;
            if (f > 0 && f < 1)
            {
                point origin((1.0-f)*origin0+f*origin1);
                point inter((1.0-f)*points[e[0]]+f*points[e[1]]);
                point originNormal
                (
                    (1.0-f)*allPointInfo[e[0]].v()
                   +f*allPointInfo[e[0]].v()
                );

                str.write(linePointRef(inter, origin));

                isoPoints.append(inter); //(inter);
                surfPoints.append(origin);      //(origin0);
                surfNormals.append(originNormal);   //(allPointInfo[e[0]].v());
            }
        }
    }


    isoPoints.shrink();
    surfPoints.shrink();
    surfNormals.shrink();


    // Slide along surface to minimise distance
    const scalar relax = 0.1;

    for (label iter = 0; iter < 10; iter++)
    {
        // Remove tangential component
        forAll(isoPoints, i)
        {
            const vector& n = surfNormals[i];
            vector d(surfPoints[i]-isoPoints[i]);
            surfPoints[i] =
                (1-relax)*surfPoints[i]
               +relax*(isoPoints[i] + (d&n)*n);
        }
        // Re-project
        List<pointIndexHit> info;
        surfMesh.findNearest
        (
            surfPoints,
            scalarField(surfPoints.size(), GREAT),
            info
        );
        vectorField normal;
        surfMesh.getNormal(info, normal);
        forAll(info, i)
        {
            surfPoints[i] = info[i].hitPoint();
            surfNormals[i] = normal[i];
        }
    }

    // Plot lines from iso back to originating surface
    {
        OBJstream str(runTime.path()/"slide_to_iso.obj");
        forAll(isoPoints, i)
        {
            str.write(linePointRef(isoPoints[i], surfPoints[i]));
        }
    }    


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
