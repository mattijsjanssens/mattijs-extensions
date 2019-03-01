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

    // // Normal topological change
    // {
    //     // Mesh changing engine.
    //     polyTopoChange meshMod(mesh);
    // 
    //     // Create mesh, return map from old to new mesh.
    //     autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);
    // 
    //     // Update fields
    //     mesh.updateMesh(map);
    // 
    //     // Optionally inflate mesh
    //     if (map().hasMotionPoints())
    //     {
    //         mesh.movePoints(map().preMotionPoints());
    //     }
    // }


    // Shuffle existing patches
    labelList oldToNew;
    label nNew;
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        // Remove zero'th
        oldToNew.setSize(pbm.size());
        forAll(pbm, patchi)
        {
            oldToNew[patchi] = patchi-1;
        }
        // Shuffle 0 to end
        oldToNew[0] = pbm.size()-1;
        // Truncate
        nNew = pbm.size()-1;
    }

DebugVar(oldToNew);
DebugVar(nNew);


    fvMeshTools::reorderPatches(mesh, oldToNew, nNew, true);

    runTime++;
    mesh.setInstance(runTime.timeName());

    Info<< "Writing mesh to " << runTime.timeName() << endl;

    mesh.write();
    


    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
