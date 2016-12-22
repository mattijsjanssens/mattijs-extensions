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
    removeFaces_agglomerate

Description
    Utility to remove faces (combines cells on both sides) using
    facePairAgglomeration.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "faceSet.H"
#include "removeFaces.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pairGAMGAgglomeration.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"


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


    while (true)
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalarField faceWeights
        (
            mag
            (
                cmptMultiply
                (
                    mesh.Sf().primitiveField()
                   /sqrt(mesh.magSf().primitiveField()),
                    vector(1, 1.01, 1.02)
                    //vector::one
                )
            )
        );

        label nCoarseCells = 0;
        labelField cellRegion
        (
            pairGAMGAgglomeration::agglomerate
            (
                nCoarseCells,
                mesh.lduAddr(),
                faceWeights
            )
        );

        label nTotalCoarseCells = returnReduce(nCoarseCells, sumOp<label>());

        if (nTotalCoarseCells == mesh.globalData().nTotalCells())
        {
            Info<< "Exiting since number of coarse cells " << nTotalCoarseCells
                << " same as the number of mesh cells "
                << mesh.globalData().nTotalCells() << nl << endl;

            break;
        }

        labelList cellRegionMaster(nCoarseCells, labelMax);
        forAll(cellRegion, celli)
        {
            label region = cellRegion[celli];
            if (cellRegionMaster[region] == labelMax)
            {
                cellRegionMaster[region] = celli;
            }
        }


        // Remove internal faces inbetween same coarse cell
        DynamicList<label> facesToRemove(mesh.nInternalFaces());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            label own = mesh.faceOwner()[facei];
            label nei = mesh.faceNeighbour()[facei];

            if (cellRegion[own] == cellRegion[nei])
            {
                facesToRemove.append(facei);
            }
        }




        // Topo changes container
        polyTopoChange meshMod(mesh);

        removeFaces faceRemover(mesh, 2);

        // Insert mesh refinement into polyTopoChange.
        faceRemover.setRefinement
        (
            facesToRemove,
            cellRegion,
            cellRegionMaster,
            meshMod
        );

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        mesh.updateMesh(morphMap);

        // Move mesh (since morphing does not do this)
        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        // Update numbering of cells/vertices.
        faceRemover.updateMesh(morphMap);


        Info<< "Coarsened from "
            << returnReduce(morphMap().nOldCells(), sumOp<label>())
            << " to "
            << nTotalCoarseCells
            << " cells" << nl << endl;

        mesh.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
