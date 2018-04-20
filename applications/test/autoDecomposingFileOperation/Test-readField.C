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

#include "unallocatedFvBoundaryMesh.H"
#include "GeometricField.H"
#include "argList.H"
#include "uVolFields.H"
#include "unallocatedFvMesh.H"
#include "unallocatedVolMesh.H"
#include "unallocatedGenericFvPatchField.H"
#include "parUnallocatedFvFieldReconstructor.H"
#include "unallocatedFvMeshTools.H"
#include "unallocatedFvFieldReconstructor.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    if (!Pstream::parRun())
    {
        // We've got IOobject so:
        //  - (undecomposed) objectRegistry (probably mesh)
        //  - (undecomposed) Time
        // Use this to load processor Times and meshes

        // Read the processor databases

DebugVar(args.path());
        label nProcs = fileHandler().nProcs(args.path(), word::null);
DebugVar(nProcs);
        PtrList<Time> databases(nProcs);
        forAll(databases, proci)
        {
            databases.set
            (
                proci,
                new Time
                (
                    Time::controlDictName,
                    args.rootPath(),
                    args.caseName()/fileName(word("processor") + name(proci))
                )
            );
        }

        PtrList<labelIOList> cellProcAddressing(nProcs);
        PtrList<labelIOList> faceProcAddressing(nProcs);
        PtrList<labelIOList> boundaryProcAddressing(nProcs);

        forAll(cellProcAddressing, proci)
        {
            cellProcAddressing.set
            (
                proci,
                new labelIOList
                (
                    IOobject
                    (
                        "cellProcAddressing",
                        runTime.constant(), //mesh.facesInstance(),
                        fvMesh::meshSubDir,
                        databases[proci],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );
        }
        forAll(faceProcAddressing, proci)
        {
            faceProcAddressing.set
            (
                proci,
                new labelIOList
                (
                    IOobject
                    (
                        "faceProcAddressing",
                        runTime.constant(), //mesh.facesInstance(),
                        fvMesh::meshSubDir,
                        databases[proci],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );
        }
        forAll(boundaryProcAddressing, proci)
        {
            boundaryProcAddressing.set
            (
                proci,
                new labelIOList
                (
                    IOobject
                    (
                        "boundaryProcAddressing",
                        runTime.constant(), //mesh.facesInstance(),
                        fvMesh::meshSubDir,
                        databases[proci],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );
        }

        // Read the (unallocated) processor meshes
        PtrList<unallocatedFvMesh> procMeshes(nProcs);

        forAll(procMeshes, proci)
        {
            Pout<< "** reading procMesh:" << proci
                << " with cellProcAddresing:"
                << cellProcAddressing[proci].size() << endl;

            procMeshes.set
            (
                proci,
                unallocatedFvMeshTools::newMesh
                (
                    IOobject
                    (
                        fvMesh::defaultRegion,
                        cellProcAddressing[proci].instance(),
                        databases[proci],
                        IOobject::MUST_READ
                    ),
                    cellProcAddressing[proci].size()
                )
            );
            Pout<< procMeshes[proci].info() << endl;
        }


        // Get mesh as unallocated
        Pout<< "** Reading baseMesh" << endl;
        autoPtr<unallocatedFvMesh> uMesh
        (
            unallocatedFvMeshTools::newMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,      // name of mesh
                    runTime.timeName(),         // start search
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        const unallocatedFvFieldReconstructor reconstructor
        (
            uMesh(),
            procMeshes,
            faceProcAddressing,
            cellProcAddressing,
            boundaryProcAddressing
        );

        // Read field on proc meshes
        PtrList<uVolScalarField> procFields(nProcs);
        forAll(procFields, proci)
        {
            const unallocatedFvMesh& procMesh = procMeshes[proci];
            procFields.set
            (
                proci,
                new uVolScalarField
                (
                    IOobject
                    (
                        "p",
                        procMesh.time().timeName(),
                        procMesh.thisDb(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    procMesh
                )
            );

            DebugVar(procFields[proci]);
        }


        // Map local field onto baseMesh
        const unallocatedFvMesh& baseMesh = uMesh();

        tmp<uVolScalarField> tfld
        (
            reconstructor.reconstructFvVolumeField
            (
                IOobject
                (
                    "p",
                    baseMesh.time().timeName(),
                    baseMesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                procFields
            )
        );

        DebugVar(tfld());
    }
    else
    {
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


        // Read procAddressing files. Deduct base mesh sizes.
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

        #include "createUnallocatedMesh.H"

        // // Read field on mesh
        // uVolScalarField baseFld
        // (
        //     IOobject
        //     (
        //         "p",
        //         baseMesh.time().timeName(),
        //         baseMesh.thisDb(),
        //         IOobject::MUST_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     baseMesh
        // );
        //uVolScalarField p
        //(
        //    IOobject
        //    (
        //        "p",
        //        baseRunTime.timeName(),
        //        baseMesh.thisDb(),
        //        IOobject::NO_READ,
        //        IOobject::AUTO_WRITE
        //    ),
        //    baseMesh,
        //    dimensionedScalar("zero", dimless, Zero),
        //    unallocatedGenericFvPatchField<scalar>::typeName
        //);

        // Mapping engine from mesh to baseMesh
        parUnallocatedFvFieldReconstructor reconstructor
        (
            baseMesh,
            mesh,
            distMap
        );

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
DebugVar(p);


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

    return 0;
}


// ************************************************************************* //
