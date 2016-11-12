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
#include "meshCutAndRemove.H"
#include "meshCutter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


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


    DynamicList<label> cutVerts;
    DynamicList<label> cutEdges;
    DynamicField<scalar> cutEdgeWeights;

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


DebugVar(cutEdges);
DebugVar(cutEdgeWeights);


    //- Construct from pattern of cuts. Detect cells to cut.
    cellCuts cuts(mesh, cutVerts, cutEdges, cutEdgeWeights);

    DebugVar(cuts.nLoops());


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.
    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.
    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);



    // Topo changes container
    polyTopoChange meshMod(mesh);

    meshCutAndRemove cutter(mesh);
    // Insert mesh refinement into polyTopoChange.
    cutter.setRefinement
    (
        0,                  //exposedPatci
        cuts,
        labelList(mesh.nCells(), 0),    //exposedPatchi
        meshMod
    );


    // Mesh change engine
    //meshCutter cutter(mesh);
    //cutter.setRefinement(cuts, meshMod);



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
    Pout<< "Writing mesh to time " << runTime.timeName() << endl;
    mesh.write();

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
