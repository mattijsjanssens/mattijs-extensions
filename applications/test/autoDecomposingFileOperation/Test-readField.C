/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Test-readField

Description
    Read volScalarField

\*---------------------------------------------------------------------------*/

#include "mapDistributeBase.H"
//#include "unallocatedFvBoundaryMesh.H"
#include "GeometricField.H"
#include "argList.H"
//#include "uSurfaceFields.H"
#include "uVolFieldsFwd.H"
//#include "uSurfaceFieldsFwd.H"
#include "unallocatedFvMesh.H"
//#include "unallocatedVolMesh.H"
//#include "unallocatedSurfaceMesh.H"
//#include "unallocatedFvsPatchField.H"
#include "unallocatedGenericFvPatchField.H"
//#include "parUnallocatedFvFieldReconstructor.H"
#include "unallocatedFvMeshTools.H"
//#include "unallocatedFvFieldReconstructor.H"
//#include "unallocatedGenericFvsPatchField.H"
//#include "uFieldReconstructor.H"
//#include "unallocatedPointMesh.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "unallocatedIOPosition.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // Load local Cloud positions
    {
        #include "createMesh.H"

        // Construct empty cloud
        IDLList<passiveParticle> dummyParticles;
        Cloud<passiveParticle> cloud(mesh, "kinematicCloud", dummyParticles);

        // Read particle positions into cloud
        unallocatedIOPosition<Cloud<passiveParticle>> ioP
        (
            IOobject
            (
                "positions2",
                cloud.time().timeName(),
                cloud,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            cloud
        );

        bool valid = ioP.headerOk();
        Istream& is = ioP.readStream("", valid);
        if (valid)
        {
            ioP.readData(is, cloud);
            ioP.close();
        }

        DebugVar(cloud.size());


        return 0;
    }


//     // Load local mesh and field
//     if (true)
//     {
//         #include "createUnallocatedMesh.H"
// 
//         uSurfaceScalarField ufld
//         (
//             IOobject
//             (
//                 "phi",
//                 mesh.time().timeName(),
//                 mesh.thisDb(),
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh
//         );
//         DebugVar(ufld);
//     }

/*
    if (!Pstream::parRun())
    {
        // Reconstruct

        #include "createUnallocatedMesh.H"

        // Read proc meshes
        const uFieldReconstructor& reconstructor = uFieldReconstructor::New
        (
            mesh.thisDb()
        );
        const PtrList<unallocatedFvMesh>& procMeshes =
            reconstructor.procMeshes();

        // Read field on proc meshes
        PtrList<uSurfaceScalarField> procFields(procMeshes.size());
        forAll(procFields, proci)
        {
            const unallocatedFvMesh& procMesh = procMeshes[proci];
            procFields.set
            (
                proci,
                new uSurfaceScalarField
                (
                    IOobject
                    (
                        "phi",
                        procMesh.time().timeName(),
                        procMesh.thisDb(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    procMesh
                )
            );
        }


        // Map local field onto baseMesh
        const unallocatedFvMesh& baseMesh = reconstructor.baseMesh();

        tmp<uSurfaceScalarField> tfld
        (
            reconstructor.reconstructor().reconstructFvSurfaceField
            (
                IOobject
                (
                    "phi",
                    baseMesh.time().timeName(),
                    baseMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                procFields
            )
        );
        DebugVar(tfld());
    }
    else
    {

const bool reconstruct = false;

        // Read procAddressing files (from runTime). Deduct base mesh sizes.
        autoPtr<mapDistributePolyMesh> distMapPtr
        (
            unallocatedFvMeshTools::readReconstructMap
            (
                IOobject
                (
                    "dummy",
                    runTime.findInstance(fvMesh::meshSubDir, "faces"),
                    fvMesh::meshSubDir,
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
        const mapDistributePolyMesh& distMap = distMapPtr();


        // Parent database
        Time baseRunTime
        (
            runTime.controlDict(),
            runTime.rootPath(),
            runTime.globalCaseName(),
            runTime.system(),
            runTime.constant(),
            false                   // enableFunctionObjects
        );
        baseRunTime.setTime(runTime);

        // Parent mesh
        autoPtr<unallocatedFvMesh> baseMeshPtr
        (
            unallocatedFvMeshTools::newMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,      // name of mesh
                    baseRunTime.timeName(),
                    baseRunTime,
                    IOobject::MUST_READ
                ),
                distMap.cellMap().constructSize()
            )
        );
        unallocatedFvMesh& baseMesh = baseMeshPtr();

        // Local mesh
        #include "createUnallocatedMesh.H"

        // Mapping engine from mesh to baseMesh
        parUnallocatedFvFieldReconstructor reconstructor
        (
            baseMesh,
            mesh,
            distMap
        );

        if (reconstruct)
        {
            // Load local field
            uVolScalarField p
            (
                IOobject
                (
                    "p",
                    mesh.time().timeName(),
                    mesh.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            );

            // Map local field onto baseMesh
            tmp<uVolScalarField> tfld
            (
                reconstructor.reconstructFvVolumeField(p)
            );

            // Write master field to parent
            tfld.ref().rename("my_p");
            {
                const bool oldParRun = Pstream::parRun();
                Pstream::parRun() = false;
                if (Pstream::master())
                {
                    DebugVar(tfld());
                    tfld().write();
                }
                Pstream::parRun() = oldParRun;
            }
        }
        else
        {
            // Decompose

            // Load base field
            uVolScalarField p
            (
                IOobject
                (
                    "p",
                    baseMesh.time().timeName(),
                    baseMesh.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                baseMesh
            );

DebugVar(p);

            // Map base field onto local mesh
            tmp<uVolScalarField> tfld
            (
                reconstructor.decomposeFvVolumeField(p)
            );
            DebugVar(tfld());
        }
    }
*/

    return 0;
}


// ************************************************************************* //
