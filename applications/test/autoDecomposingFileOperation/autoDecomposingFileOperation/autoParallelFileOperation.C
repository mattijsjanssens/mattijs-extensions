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
#include "addToRunTimeSelectionTable.H"
#include "mapDistributePolyMesh.H"
#include "unallocatedFvMesh.H"
#include "unallocatedFvMeshTools.H"
#include "parUnallocatedFvFieldReconstructor.H"
#include "uVolFields.H"
#include "unallocatedFvMeshObject.H"
#include "streamReconstructor.H"
#include "unthreadedInitialise.H"

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

    // Mark as not needing threaded mpi
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        unthreadedInitialise,
        word,
        autoParallel
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
        // Install default file handler
        storeFileHandler defaultOp;

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

        // Use whatever we read
        const_cast<Time&>(runTime).setTime(baseRunTimePtr_());
    }
//    else
//    {
//        baseRunTimePtr_().setTime(runTime);
//    }
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
        // Install default file handler
        storeFileHandler defaultOp;

        procFacesInstance_ = runTime.findInstance(fvMesh::meshSubDir, "faces");

        Info<< "Reading decomposition addressing from " << procFacesInstance_
            << nl << endl;
        distMapPtr_ = unallocatedFvMeshTools::readReconstructMap
        (
            IOobject
            (
                "dummy",
                procFacesInstance_,
                fvMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );


        const Time& baseTime = baseRunTime(runTime);

        Info<< "Reading base mesh for time = "
            << baseTime.timeName() << nl << endl;
        baseMeshPtr_ = unallocatedFvMeshTools::newMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,      // name of mesh
                baseTime.timeName(),
                baseTime,
                IOobject::MUST_READ
            ),
            distMapPtr_().cellMap().constructSize()
        );
    }
    return baseMeshPtr_();
}


const Foam::unallocatedFvMesh&
Foam::fileOperations::autoParallelFileOperation::mesh(const Time& runTime) const
{
    if (!meshPtr_.valid())
    {
        // Install default file handler
        storeFileHandler defaultOp;

        Info<< "Reading current mesh" << nl << endl;
        meshPtr_ = Foam::unallocatedFvMeshTools::newMesh
        (
            Foam::IOobject
            (
                Foam::fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        );
    }
    return meshPtr_();
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
        Pout<< indent
            << "autoParallelFileOperation::filePath :"
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
        Pout<< indent
            << "autoParallelFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::instantList Foam::fileOperations::autoParallelFileOperation::findTimes
(
    const fileName& dir,
    const word& constantName
) const
{
    instantList times;
    if (!Pstream::parRun())
    {
        times = uncollatedFileOperation::findTimes(dir, constantName);
    }
    else
    {
        if (debug)
        {
            Pout<< indent
                << "autoParallelFileOperation::findTimes :"
                << " Searching in parent " << dir/".."
                << endl;
        }
        times = uncollatedFileOperation::findTimes(dir/"..", constantName);
    }

    if (debug)
    {
        Pout<< indent
            << "autoParallelFileOperation::findTimes :"
            << " Returning from time searching:" << endl << indent
            << "    dir   :" << dir << endl << indent
            << "    times :";
        if (times.size() >= 2)
        {
            Pout<< times[0].name() << " .. " << times.last().name();
        }
        else
        {
            Pout<< times;
        }
        Pout<< endl << endl;
    }
    return times;
}


void Foam::fileOperations::autoParallelFileOperation::setTime
(
    const Time& runTime
) const
{
    if (Pstream::parRun() && runTime.processorCase())
    {
        if (debug)
        {
            Pout<< indent
                << "autoParallelFileOperation::setTime :"
                << " Setting base time to " << runTime.timeName()
                << endl;
        }
        const Time& baseTime = baseRunTime(runTime);
        const_cast<Time&>(baseTime).setTime(runTime);
    }
}


Foam::fileNameList
Foam::fileOperations::autoParallelFileOperation::readObjects
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    word& newInstance
) const
{
    fileNameList objects;

    if (!Pstream::parRun())
    {
        objects = uncollatedFileOperation::readObjects
        (
            db,
            instance,
            local,
            newInstance
        );
    }
    else
    {
        // tbd: lagrangian. dbDir="", local = "lagrangian/KinematicCloud"
        fileName path
        (
            filePath
            (
                db.rootPath()
               /db.time().globalCaseName()
               /instance
               /db.dbDir()
               /local
            )
        );

        if (Foam::isDir(path))
        {
            newInstance = instance;
            objects = Foam::readDir(path, fileName::FILE);

            if (debug)
            {
                Pout<< indent
                    << "autoParallelFileOperation::readObjects :"
                    << " Returning parent directory searching:"
                    << endl << indent
                    << "    path     :" << path << endl << indent
                    << "    objects  :" << objects << endl << endl;
            }

            return objects;
        }
        else
        {
            objects = uncollatedFileOperation::readObjects
            (
                db,
                instance,
                local,
                newInstance
            );
        }
    }

    if (debug)
    {
        Pout<< indent
            << "autoParallelFileOperation::readObjects :"
            << " Returning from directory searching:" << endl << indent
            << "    path     :" << db.path(instance, db.dbDir()/local)
            << endl << indent
            << "    objects  :" << objects << endl << indent
            << "    newInst  :" << newInstance << endl << endl;
    }
    return objects;
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
                << " Searching for handler for type:" << type
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
                Pout<< indent
                    << "io.objectPath   :" << io.objectPath() << nl
                    << indent
                    << "filePath        :" << objPath << nl
                    << indent
                    << "parentObjectPath:" << parentObjectPath << endl;
            }

            if (io.objectPath() != objPath && objPath == parentObjectPath)
            {
                const Time& runTime = io.time();

                // Install default file handler
                storeFileHandler defaultOp;

                Pout<< incrIndent;

                // Read local mesh
                const unallocatedFvMesh& procMesh = mesh(runTime);
                // Read undecomposed mesh. Read procAddressing files
                // (from runTime).
                const unallocatedFvMesh& baseUMesh = baseMesh(runTime);
                const mapDistributePolyMesh& distMap = distMapPtr_();

                // Mapping engine from mesh to baseMesh
                const parUnallocatedFvFieldReconstructor reconstructor
                (
                    baseUMesh,
                    procMesh,
                    distMap
                );

                IOobject baseIO
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    baseUMesh.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                );

                OStringStream os(IOstream::BINARY);

                if (debug)
                {
                    Pout<< "autoParallelFileOperation::read :"
                        << " decompose and writing:" << baseIO.objectPath()
                        << endl;
                }
                bool ok = typeReconstructor().decompose
                (
                    reconstructor,
                    baseUMesh,
                    baseIO,
                    procMesh,
                    io,
                    false,              // no face flips. Tbd.
                    os
                );

                Pout<< decrIndent;

                if (ok)
                {
                    IStringStream is(os.str(), IOstream::BINARY);

                    // Read field from stream
                    ok = io.readData(is);
                    io.close();

                    if (debug)
                    {
                        Pout<< indent
                            << "autoParallelFileOperation::read :"
                            << " sucessfully decomposed " << io.objectPath()
                            << endl;
                    }
                }
                else
                {
                    if (debug)
                    {
                        Pout<< indent
                            << "autoParallelFileOperation::read :"
                            << " ** failed decomposing " << io.objectPath()
                            << endl;
                    }
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
        if (debug)
        {
            Pout<< indent
                << "autoParallelFileOperation::writeObject :"
                << " Searching for handler for type:" << io.type()
                << " of object: " << io.objectPath() << endl;
        }

        autoPtr<streamReconstructor> typeReconstructor
        (
            streamReconstructor::New(io.type())
        );

        if (typeReconstructor.valid())
        {
            const Time& runTime = io.time();

            Pout<< incrIndent;

            // Install default file handler
            storeFileHandler defaultOp;

            // Read local mesh
            const unallocatedFvMesh& procMesh = mesh(runTime);
            // Read undecomposed mesh. Read procAddressing files
            // (from runTime).
            const unallocatedFvMesh& baseUMesh = baseMesh(runTime);
            const mapDistributePolyMesh& distMap = distMapPtr_();

            // Mapping engine from mesh to baseMesh
            const parUnallocatedFvFieldReconstructor reconstructor
            (
                baseUMesh,
                procMesh,
                distMap
            );

            if (debug)
            {
                Pout<< indent
                    << "autoParallelFileOperation::writeObject :"
                    << " reconstructing and writing:" << io.name()
                    << " with handler for type:" << io.type()
                    << endl;
            }

            Pout<< decrIndent;

            typeReconstructor().reconstruct
            (
                reconstructor,
                io,
                false,          // no face flips. Tbd.
                fmt,
                ver,
                cmp
            );
        }
        else if
        (
            io.instance() == io.time().timeName()
         && io.local() == "uniform"
        )
        {
            //Pout<< indent << "** uniform:" << io.objectPath() << endl;
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
