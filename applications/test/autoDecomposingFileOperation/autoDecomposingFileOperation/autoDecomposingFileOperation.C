/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "autoDecomposingFileOperation.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvFieldDecomposer.H"
#include "addToRunTimeSelectionTable.H"
#include "processorMeshes.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(autoDecomposingFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        autoDecomposingFileOperation,
        word
    );


    class installFileOp
    {
        word oldHandler_;

    public:

        installFileOp()
        {
            // Install autoDecomposing as fileHandler
            oldHandler_ = fileHandler().type();
            autoPtr<fileOperation> handler
            (
                new autoDecomposingFileOperation(true)
            );
            Foam::fileHandler(handler);
        }

        ~installFileOp()
        {
            //// Restore old fileHandler
            //autoPtr<fileOperation> handler
            //(
            //    fileOperation::New
            //    (
            //        oldHandler_,
            //        true
            //    )
            //);
            autoPtr<fileOperation> handler(nullptr);
            Foam::fileHandler(handler);
        }
    };
    installFileOp installFileOp_;
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fvMesh& Foam::fileOperations::autoDecomposingFileOperation::baseMesh
(
    const Time& runTime
) const
{
    if (!baseMeshPtr_.valid())
    {
        Info<< "Create base time\n" << endl;
        baseRunTimePtr_.reset
        (
            new Time
            (
                runTime.controlDict(),
                runTime.rootPath(),
                runTime.globalCaseName(),
                runTime.system(),
                runTime.constant(),
                false                   // enableFunctionObjects
            )
        );

        Info<< "Create base mesh for time = "
            << runTime.timeName() << nl << endl;
        baseMeshPtr_.reset
        (
            new fvMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    baseRunTimePtr_().timeName(),
                    baseRunTimePtr_(),
                    IOobject::MUST_READ
                )
            )
        );
    }
    return baseMeshPtr_();
}


const Foam::fvFieldDecomposer&
Foam::fileOperations::autoDecomposingFileOperation::decomposer
(
    const IOobject& io
) const
{
    // Force reading of the decomposition maps
    if (!decomposerPtr_.valid())
    {
        const fvMesh& procMesh = dynamic_cast<const fvMesh&>(io.db());
        // Force reading of the baseMesh
        const fvMesh& undecomposedMesh = baseMesh(io.time());

        if (debug)
        {
            Pout<< "Loading maps\n" << endl;
        }
        cellAddressingPtr_.reset
        (
            new labelIOList
            (
                IOobject
                (
                    "cellProcAddressing",
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
        faceAddressingPtr_.reset
        (
            new labelIOList
            (
                IOobject
                (
                    "faceProcAddressing",
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
        boundaryAddressingPtr_.reset
        (
            new labelIOList
            (
                IOobject
                (
                    "boundaryProcAddressing",
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
        if (debug)
        {
            Pout<< "Constructing decomposer\n" << endl;
        }
        decomposerPtr_.reset
        (
            new fvFieldDecomposer
            (
                undecomposedMesh,
                procMesh,
                faceAddressingPtr_(),
                cellAddressingPtr_(),
                boundaryAddressingPtr_()
            )
        );
    }
    return decomposerPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::autoDecomposingFileOperation::
autoDecomposingFileOperation
(
    const bool verbose
)
:
    uncollatedFileOperation(false)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::autoDecomposingFileOperation::
~autoDecomposingFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::autoDecomposingFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io
) const
{
    if (debug)
    {
        Pout<< "autoDecomposingFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    // Try uncollated searching
    fileName objPath = uncollatedFileOperation::filePath(checkGlobal, io);

    // If not found and parallel check parent
    if (objPath.empty() && checkGlobal && io.time().processorCase())
    {
        fileName parentObjectPath =
            io.rootPath()/io.time().globalCaseName()
           /io.instance()/io.db().dbDir()/io.local()/io.name();

        if (isFile(parentObjectPath))
        {
            objPath = parentObjectPath;
        }
    }

    if (debug)
    {
        Pout<< "autoDecomposingFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


bool Foam::fileOperations::autoDecomposingFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& type
) const
{
    bool ok = true;

    if
    (
        Pstream::parRun()
     && (
            type == volScalarField::typeName
         || type == volVectorField::typeName
         || type == volSphericalTensorField::typeName
         || type == volSymmTensorField::typeName
         || type == volTensorField::typeName
        )
    )
    {
        // Set flag for e.g. codeStream
        const bool oldGlobal = io.globalObject();
        io.globalObject() = masterOnly;
        // If codeStream originates from dictionary which is
        // not IOdictionary we have a problem so use global
        //const bool oldFlag = regIOobject::masterOnlyReading;
        //regIOobject::masterOnlyReading = masterOnly;


        // Find file, check in parent directory
        fileName objPath = filePath(true, io);

        // Check if the file comes from the parent path
        fileName parentObjectPath =
            io.rootPath()/io.time().globalCaseName()
           /io.instance()/io.db().dbDir()/io.local()/io.name();

        if (debug)
        {
            Pout<< "io.objectPath   :" << io.objectPath() << endl;
            Pout<< "filePath        :" << objPath << endl;
            Pout<< "parentObjectPath:" << parentObjectPath << endl;
        }

        if (io.objectPath() != objPath && objPath == parentObjectPath)
        {
            // Force reading of the baseMesh
            const fvMesh& undecomposedMesh = baseMesh(io.time());

            IOobject parentIO
            (
                io.name(),
                io.instance(),
                io.local(),
                undecomposedMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            OStringStream os;

            decomposeAndWrite<volScalarField>(io, parentIO, type, os);
            decomposeAndWrite<volVectorField>(io, parentIO, type, os);
            decomposeAndWrite<volSphericalTensorField>(io, parentIO, type, os);
            decomposeAndWrite<volSymmTensorField>(io, parentIO, type, os);
            decomposeAndWrite<volTensorField>(io, parentIO, type, os);

            IStringStream is(os.str());

            // Read field from stream
            ok = io.readData(is);
            io.close();
        }
        else
        {
            ok = io.readData(io.readStream(type));
            io.close();
        }

        // Restore flags
        io.globalObject() = oldGlobal;
        //regIOobject::masterOnlyReading = oldFlag;
    }
    else
    {
        ok = uncollatedFileOperation::read
        (
            io,
            masterOnly,
            format,
            type
        );
    }

    return ok;
}


bool Foam::fileOperations::autoDecomposingFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    bool ok = true;

    if (!valid)
    {
        return ok;
    }

    if
    (
        Pstream::parRun()
     && (
            io.type() == volScalarField::typeName
         //|| io.type() == volVectorField::typeName
         //|| io.type() == volSphericalTensorField::typeName
         //|| io.type() == volSymmTensorField::typeName
         //|| io.type() == volTensorField::typeName
        )
    )
    {
        if (Pstream::master())
        {
            bool oldParRun = UPstream::parRun();
            UPstream::parRun() = false;

            // Load undecomposed mesh
            const fvMesh& base = baseMesh(io.time());

            // Create the processor databases
            PtrList<Time> databases(Pstream::nProcs());

            forAll(databases, proci)
            {
                databases.set
                (
                    proci,
                    new Time
                    (
                        io.time().controlDict(),
                        io.time().rootPath(),
                        io.time().globalCaseName()
                       /fileName(word("processor") + name(proci)),
                        io.time().system(),
                        io.time().constant(),
                        false                   // enableFunctionObjects
                    )
                );
            }

            // Read all meshes and addressing to reconstructed mesh
            processorMeshes procMeshes(databases, io.db().name());

            UPstream::parRun() = oldParRun;
        }

        // Get the fields locally
        PstreamBuffers pBufs(Pstream::nonBlocking);
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOpstream os(Pstream::masterNo());
            io.writeData(os);
        }
        pBufs.finishedSends();

        PtrList<volScalarField> procFields(Pstream::nProcs());

        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            const fvMesh& procMesh = procMeshes.meshes()[proci];

            IOobject procIO(io, procMesh);

            UIPstream is(Pstream::masterNo());
            procFields.set
            (
                proci,
                new volScalarField
                (
                    procIO,
                    procMesh,
                    dictionary(is)
                )
            );
        }

        tmp<volScalarField> tfld
        (
            fvReconstructor.reconstructFvVolumeField
            (
                baseIO,
                procFields
            )
        );
        tfld().writeObject(fmt, ver, cmp);
    }
//XXXX



    }
    else
    {
        ok = uncollatedFileOperation::writeObject
        (
            io,
            fmt,
            ver,
            cmp,
            valid
        );
    }
    return ok;
}


// ************************************************************************* //
