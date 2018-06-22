/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "autoParallelFileOperation.H"
#include "Time.H"
#include "fvMesh.H"
//#include "fvFieldDecomposer.H"
#include "addToRunTimeSelectionTable.H"
//#include "processorMeshes.H"
//#include "fvFieldReconstructor.H"

#include "mapDistributePolyMesh.H"
#include "unallocatedFvMesh.H"
#include "unallocatedFvMeshTools.H"
#include "parUnallocatedFvFieldReconstructor.H"
#include "uVolFields.H"
#include "unallocatedFvMeshObject.H"

#include "streamReconstructor.H"
#include "parUnallocatedFvFieldReconstructor.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(autoParallelFileOperation, 0);
    addRemovableToRunTimeSelectionTable
    (
        fileOperation,
        autoParallelFileOperation,
        word
    );

    class installFileOp
    {
    public:

        installFileOp()
        {
            // Install autoDecomposing as fileHandler
            //autoPtr<fileOperation> handler
            //(
            //    new autoParallelFileOperation(true)
            //);
            //Foam::fileHandler(handler);
        }

        ~installFileOp()
        {
            if (fileHandler().type() == autoParallelFileOperation::typeName)
            {
                autoPtr<fileOperation> handler(nullptr);
                Foam::fileHandler(handler);
            }
        }
    };
    installFileOp installAutoParallel_;
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Time& Foam::fileOperations::autoParallelFileOperation::baseRunTime
(
    const Time& runTime
) const
{
    if (!baseRunTimePtr_.valid())
    {
        Info<< "Creating base time\n" << endl;
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
    }
    baseRunTimePtr_().setTime(runTime);
    return baseRunTimePtr_();
}


const Foam::unallocatedFvMesh&
Foam::fileOperations::autoParallelFileOperation::baseMesh
(
    const Time& runTime
) const
{
    if (!baseMeshPtr_.valid())
    {
        Info<< "Reading decomposition addressing" << nl << endl;
        distMapPtr_ = unallocatedFvMeshTools::readReconstructMap
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
        );

        Info<< "Creating base mesh for time = "
            << runTime.timeName() << nl << endl;
        baseMeshPtr_ = unallocatedFvMeshTools::newMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,      // name of mesh
                baseRunTime(runTime).timeName(),
                baseRunTime(runTime),
                IOobject::MUST_READ
            ),
            distMapPtr_().cellMap().constructSize()
        );
    }
    return baseMeshPtr_();
}


// const Foam::parUnallocatedFvFieldReconstructor&
// Foam::fileOperations::autoParallelFileOperation::reconstructor
// (
//     const unallocatedFvMesh& mesh
// ) const
// {
//     if (!reconstructorPtr_.valid())
//     {
//         Info<< "Creating reconstructor" << nl << endl;
//         reconstructorPtr_ = new parUnallocatedFvFieldReconstructor
//         (
//             baseMesh(),
//             mesh,
//             distMapPtr_()
//         );
//     }
//     return reconstructorPtr_();
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::autoParallelFileOperation::
autoParallelFileOperation
(
    const bool verbose
)
:
    uncollatedFileOperation(false)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName << nl
            << "         (reconstructs on-the-fly in parallel operation)"
            << endl;
    }

    if (regIOobject::fileModificationChecking == regIOobject::timeStampMaster)
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
    }
    else if
    (
        regIOobject::fileModificationChecking
     == regIOobject::inotifyMaster
    )
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify"
                << endl;
        }
        regIOobject::fileModificationChecking = regIOobject::inotify;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::autoParallelFileOperation::
~autoParallelFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::autoParallelFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& typeName
) const
{
    if (debug)
    {
        Pout<< "autoParallelFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    // Try uncollated searching
    fileName objPath = uncollatedFileOperation::filePath
    (
        checkGlobal,
        io,
        typeName
    );

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
        Pout<< "autoParallelFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


bool Foam::fileOperations::autoParallelFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& type
) const
{
    bool ok = true;

    if (Pstream::parRun())
    {
        if (debug)
        {
            Pout<< indent
                << "autoParallelFileOperation::read :"
                << " Searching for reconstructor for type:" << type
                << " of object: " << io.objectPath() << endl;
        }
        autoPtr<streamReconstructor> typeReconstructor
        (
            streamReconstructor::New(type)
        );

        if (typeReconstructor.valid())
        {
            // Set flag for e.g. codeStream
            const bool oldGlobal = io.globalObject();
            io.globalObject() = masterOnly;
            // If codeStream originates from dictionary which is
            // not IOdictionary we have a problem so use global
            //const bool oldFlag = regIOobject::masterOnlyReading;
            //regIOobject::masterOnlyReading = masterOnly;


            // Find file, check in parent directory
            fileName objPath = filePath(true, io, type);

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
                const Time& runTime = io.time();

                // Read procAddressing files (from runTime). Deduct base
                // mesh sizes.
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
                const parUnallocatedFvFieldReconstructor reconstructor
                (
                    baseMesh,
                    mesh,
                    distMap
                );

                IOobject baseIO
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    baseMesh.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                );

                OStringStream os(IOstream::BINARY);

                if
                (
                    typeReconstructor().decompose
                    (
                        reconstructor,
                        baseMesh,
                        baseIO,
                        mesh,
                        io,
                        os
                    )
                )
                {
                    IStringStream is(os.str(), IOstream::BINARY);

                    // Read field from stream
                    ok = io.readData(is);
                    io.close();
                }
                else
                {
                    return false;
                }
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
            ok = io.readData(io.readStream(type));
            io.close();
        }
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


bool Foam::fileOperations::autoParallelFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    bool ok = true;

    if (Pstream::parRun())
    {
//XXX
/*
        autoPtr<streamReconstructor> typeReconstructor
        (
            streamReconstructor::New(type)
        );

        if (typeReconstructor.valid())
        {
            typeReconstructor.reconstruct
            (
                reconstructor,
                io,
                fmt,
                ver,
                cmp
            );
*/
//XXX


        if
        (
            io.type() == volScalarField::typeName
         || io.type() == volVectorField::typeName
         || io.type() == volSphericalTensorField::typeName
         || io.type() == volSymmTensorField::typeName
         || io.type() == volTensorField::typeName
         || io.type() == surfaceScalarField::typeName
        )
        {
            if (!valid)
            {
                return ok;
            }

            if (io.type() == volScalarField::typeName)
            {
                return reconstructAndWrite<volScalarField>(io, fmt, ver, cmp);
            }
            else if (io.type() == volVectorField::typeName)
            {
                return reconstructAndWrite<volVectorField>(io, fmt, ver, cmp);
            }
            else if (io.type() == volSphericalTensorField::typeName)
            {
                return reconstructAndWrite<volSphericalTensorField>
                (
                    io,
                    fmt,
                    ver,
                    cmp
                );
            }
            else if (io.type() == volSymmTensorField::typeName)
            {
                return reconstructAndWrite<volSymmTensorField>
                (
                    io,
                    fmt,
                    ver,
                    cmp
                );
            }
            else if (io.type() == volTensorField::typeName)
            {
                return reconstructAndWrite<volTensorField>(io, fmt, ver, cmp);
            }
            else if (io.type() == surfaceScalarField::typeName)
            {
                return reconstructAndWrite2<surfaceScalarField>
                (
                    io,
                    fmt,
                    ver,
                    cmp
                );
            }
            else
            {
                FatalErrorInFunction << "Problem" << exit(FatalError);
                return false;
            }
        }
        else if
        (
            io.instance() == io.time().timeName()
         && io.local() == "uniform"
        )
        {
            Pout<< "** uniform:" << io.objectPath() << endl;
            if (valid)
            {
                // Copy of fileOperation::writeObject but with parent path

                fileName pathName
                (
                    io.rootPath()
                   /io.time().globalCaseName()
                   /io.instance()
                   /io.db().dbDir()
                   /io.local()
                   /io.name()
                );

                mkDir(pathName.path());

                autoPtr<Ostream> osPtr
                (
                    NewOFstream
                    (
                        pathName,
                        fmt,
                        ver,
                        cmp
                    )
                );

                if (!osPtr.valid())
                {
                    return false;
                }

                Ostream& os = osPtr();

                // If any of these fail, return (leave error handling to
                // Ostream class)
                if (!os.good())
                {
                    return false;
                }

                if (!io.writeHeader(os))
                {
                    return false;
                }

                // Write the data to the Ostream
                if (!io.writeData(os))
                {
                    return false;
                }

                IOobject::writeEndDivider(os);
            }
            return true;
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
    }
    else
    {
Pout<< "** non-par:" << io.objectPath() << endl;

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
