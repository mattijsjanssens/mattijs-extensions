/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Test-pointMesh

Description

\*---------------------------------------------------------------------------*/

#include "MeshObject.H"
#include "fvMeshTools.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "volMesh.H"
#include "fvMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "IOobjectList.H"
#include "polyTopoChange.H"
#include "wallPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"
 
    const bool fields = true;

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    #include "readVolFields.H"
    #include "readSurfaceFields.H"
    #include "readPointFields.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Delete existing patch
    {
        const label deletePatchi = 0;

        Info<< "Removing patch " << deletePatchi
            << " name " << pbm[deletePatchi].name() << endl;

        // Remove deletePatchi
        labelList newToOld;
        newToOld.setSize(pbm.size()-1);
        label newi = 0;
        forAll(pbm, patchi)
        {
            if (patchi == deletePatchi)
            {
            }
            else
            {
                newToOld[newi++] = patchi;
            }
        }
        mesh.reorderPatches(newToOld, true);

        runTime++;
        mesh.setInstance(runTime.timeName());
        Info<< "Writing removed-patch mesh to " << runTime.timeName() << endl;
        mesh.write();
    }



    // Insert new (global) patch
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        label insertPatchi = pbm.size();
        forAllReverse(pbm, patchi)
        {
            if (!isA<processorPolyPatch>(pbm[patchi]))
            {
                insertPatchi = patchi+1;
                break;
            }
        }

        wallPolyPatch pp
        (
            "myPatch",
            0,              // dummy
            0,              // dummy
            0,              // dummy
            mesh.boundaryMesh(),
            wallPolyPatch::typeName
        );

        Info<< "Inserting patch " << pp.name()
            << " at location " << insertPatchi << endl;

        dictionary patchFieldDict;

        mesh.addPatch
        (
            insertPatchi,
            pp,
            patchFieldDict,
            fvPatchField<scalar>::calculatedType(),   //defaultPatchFieldType,
            true                // validBoundary
        );

        runTime++;
        mesh.setInstance(runTime.timeName());
        Info<< "Writing inserted patch mesh to " << runTime.timeName() << endl;
        mesh.write();
    }




    // Append local patch
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        label insertPatchi = pbm.size();

        wallPolyPatch pp
        (
            "localPatch",
            0,              // dummy
            0,              // dummy
            0,              // dummy
            mesh.boundaryMesh(),
            wallPolyPatch::typeName
        );

        Info<< "Inserting patch " << pp.name()
            << " at location " << insertPatchi << endl;

        dictionary patchFieldDict;

        mesh.addPatch
        (
            insertPatchi,
            pp,
            patchFieldDict,
            fvPatchField<scalar>::calculatedType(),   //defaultPatchFieldType,
            false                // validBoundary
        );

        runTime++;
        mesh.setInstance(runTime.timeName());
        Info<< "Writing LOCAL inserted patch mesh to "
            << runTime.timeName() << endl;
        mesh.write();
    }



    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
