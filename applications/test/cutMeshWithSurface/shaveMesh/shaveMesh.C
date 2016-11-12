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


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
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
#include "OBJstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addOverwriteOption.H"
    argList::validArgs.append("patch");

    #include "setRootCase.H"
    #include "createTime.H"


    runTime.functionObjects().off();
    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs[0], 0);

    #include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const fileName patchName = args[1];
    const bool overwrite = args.optionFound("overwrite");

    label patchi = mesh.boundaryMesh().findPatchID(patchName);

    DebugVar(patchi);

    const polyPatch& pp = mesh.boundaryMesh()[patchi];

    const labelList& own = mesh.faceOwner();
    const labelList& meshPoints = pp.meshPoints();
    const labelListList& pointFaces = pp.pointFaces();
    const labelListList& pointEdges = pp.pointEdges();
    const edgeList& edges = pp.edges();



    DynamicList<label> cutCells;
    DynamicList<label> cutVerts;
    DynamicList<label> cutEdges;
    DynamicField<scalar> cutEdgeWeights;

    {
        OBJstream str(runTime.path()/"cutVerts.obj");

        PackedBoolList isCellCut(mesh.nCells());
        PackedBoolList isVertCut(mesh.nCells());


        forAll(pointFaces, ppi)
        {
            const labelList& pFaces = pointFaces[ppi];
            if (pFaces.size() == 3)
            {
                // Check that the three faces are all on the same cell
                label celli = own[pp.start()+pFaces[0]];
                if
                (
                    own[pp.start()+pFaces[1]] == celli
                 && own[pp.start()+pFaces[2]] == celli
                )
                {
                    label pointi = meshPoints[ppi];
                    //DebugVar(mesh.points()[pointi]);

                    if (isCellCut.set(celli))
                    {
                        cutCells.append(celli);
                    }
                    const labelList& pEdges = pointEdges[ppi];
                    forAll(pEdges, i)
                    {
                        label edgei = pEdges[i];
                        const edge& e = edges[edgei];
                        label otherPointi = e.otherVertex(ppi);

                        if (isVertCut.set(meshPoints[otherPointi]))
                        {
                            cutVerts.append(meshPoints[otherPointi]);
                        }

                        str.write
                        (
                            linePointRef
                            (
                                mesh.points()[pointi],
                                mesh.points()[meshPoints[otherPointi]]
                            )
                        );
                    }
                }
            }
        }
    }



    //- Construct from pattern of cuts. Detect cells to cut.
    cellCuts cuts(mesh, cutCells, cutVerts, cutEdges, cutEdgeWeights);

    DebugVar(cuts.nLoops());

    {
        OBJstream str(runTime.path()/"anchors.obj");
        const labelListList& cellLoops = cuts.cellLoops();
        const labelListList& cellAnchorPoints = cuts.cellAnchorPoints();

        forAll(cellLoops, celli)
        {
            if (cellLoops[celli].size())
            {
                const labelList& anchors = cellAnchorPoints[celli];

                if (anchors.size() == 1)
                {
                    cuts.flip(celli);
                }

                const point& cc = mesh.cellCentres()[celli];
                forAll(anchors, i)
                {
                    const point& anchorPt = mesh.points()[anchors[i]];
                    str.write(linePointRef(cc, anchorPt));
                }
            }
        }
    }



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
        labelList(mesh.nCells(), patchi),    //exposedPatchi
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
