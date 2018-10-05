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

\*---------------------------------------------------------------------------*/

#include "autoReconstructFileOperation.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "unthreadedInitialise.H"
#include "streamReconstructor.H"
#include "cloud.H"
#include "dummyISstream.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(autoReconstructFileOperation, 0);
    addRemovableToRunTimeSelectionTable
    (
        fileOperation,
        autoReconstructFileOperation,
        word
    );

    // Mark as not needing threaded mpi
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        unthreadedInitialise,
        word,
        autoReconstruct
    );


    class autoReconstructInstall
    {
    public:

        autoReconstructInstall()
        {}

        ~autoReconstructInstall()
        {
            // Uninstall me as a file handler
            if
            (
                fileOperation::fileHandlerPtr_.valid()
             && (
                    fileOperation::fileHandlerPtr_().type()
                 == autoReconstructFileOperation::typeName
                )
            )
            {
                // Install 'simple' file handler to avoid e.g. communicator
                // allocation
                autoPtr<fileOperation> handler
                (
                    new uncollatedFileOperation(false)
                );
                Foam::fileHandler(handler);
            }
        }
    };
    autoReconstructInstall autoReconstructInstallObject_;


    //- Less function class used in sorting instants
    class lessOrConstant
    {
        const UList<instant>& times_;

    public:

        lessOrConstant(const UList<instant>& times)
        :
            times_(times)
        {}

        bool operator()(const label a, const label b) const
        {
            const instant& ia = times_[a];
            const instant& ib = times_[b];

            if (ia.name() == "constant")
            {
                return true;
            }
            else if (ib.name() == "constant")
            {
                return false;
            }
            else
            {
                return ia.value() < ib.value();
            }
        }

        bool equal(const label a, const label b) const
        {
            const instant& ia = times_[a];
            const instant& ib = times_[b];

            if (ia.name() == "constant")
            {
                return (ib.name() == "constant");
            }
            else if (ib.name() == "constant")
            {
                return false;
            }
            else
            {
                return ia.value() == ib.value();
            }
        }
    };

}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileOperations::autoReconstructFileOperation::haveProcPath
(
    const IOobject& io,
    fileName& procObjectPath
) const
{
    if (Pstream::parRun())
    {
        return false;
    }
    else
    {
        // tbd. lagrangian: scan all processor directories.
        //label numProcs = fileOperation::nProcs(io.time().path(), word::null);
        label numProcs = 1;

        for (label proci = 0; proci < numProcs; proci++)
        {
            procObjectPath = filePath
            (
                io.rootPath()
               /io.caseName()
               /"processor" + Foam::name(proci)
               /io.instance()
               /io.db().dbDir()
               /io.local()
               /io.name()
            );
            if (procObjectPath.size())
            {
                return true;
            }
        }
        return false;
    }
}


Foam::fileName
Foam::fileOperations::autoReconstructFileOperation::equivalentLagrangian
(
    const fileName& dir
) const
{
    std::string::size_type pos = dir.find(cloud::prefix);
    if (pos == string::npos)
    {
        return dir;
    }
    else
    {
        const wordList components(dir.components());
        label index = findIndex(components, cloud::prefix);

        if (index == -1 || index < 1)
        {
            return dir;
        }
        else
        {
            // Reconstruct path
            fileName procDir;
            if (dir[0] == '/')
            {
                procDir = '/';
            }
            forAll(components, i)
            {
                if (i != 0)
                {
                    procDir += '/';
                }
                if (i == index-1)
                {
                    procDir += "processor0/";
                }
                procDir += components[i];
            }

            //Pout<< "was:" << dir << " now:" << procDir << endl;

            return procDir;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::autoReconstructFileOperation::
autoReconstructFileOperation
(
    const bool verbose
)
:
    uncollatedFileOperation(false)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName << nl
            << "         (reconstructs on-the-fly in serial operation)" << endl;
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

    // Construct basic file handler
    basicFileHandler_ = fileOperation::New
    (
        uncollatedFileOperation::typeName,
        true
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::autoReconstructFileOperation::
~autoReconstructFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::autoReconstructFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& typeName
) const
{
    if (debug)
    {
        Pout<< indent
            << "autoReconstructFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " typeName:" << typeName
            << " checkGlobal:" << checkGlobal << endl;
    }

    // Here we want to make sure we detect any Field files in processor
    // directories but all other files (dictionaries, mesh files) in the
    // normal location.
    // The problem is that e.g. IOobjectList does a readDir followed
    // by a typeHeaderOk checking with 'labelIOList' as typeName. So
    // here we cannot checked for the wanted type (eg. volScalarField). Instead
    // we cache the object names from a previous call to readObjects and
    // check if it is one of them.

     // Cached field file?
     HashPtrTable<fileNameList, fileName>::const_iterator tmFnd =
         procObjects_.find(io.path());

    if
    (
        !checkGlobal
      && tmFnd != procObjects_.end()
      && findIndex(*tmFnd(), io.name()) != -1
    )
    {
        if (debug)
        {
            Pout<< indent
                << "autoReconstructFileOperation::filePath :"
                << " Found cached object" << endl << indent
                << "    name      :" << io.name() << endl << indent
                << "    at path   :" << io.path() << endl << endl;
        }

        fileName procPath;
        if (haveProcPath(io, procPath))
        {
            if (debug)
            {
                Pout<< indent
                    << "autoReconstructFileOperation::filePath :"
                    << " Returning from processor searching:" << endl << indent
                    << "    objectPath:" << io.objectPath()
                    << " type:" << typeName
                    << endl << indent
                    << "    filePath  :" << procPath << endl << endl;
            }
            return procPath;
        }
    }

    // Try uncollated searching
    const fileName objPath
    (
        uncollatedFileOperation::filePath
        (
            checkGlobal,
            io,
            typeName
        )
    );

    if (debug)
    {
        Pout<< indent
            << "autoReconstructFileOperation::filePath :"
            << " Returning from file searching:" << endl << indent
            << "    objectPath:" << io.objectPath() << " type:" << typeName
            << endl << indent
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileNameList Foam::fileOperations::autoReconstructFileOperation::readDir
(
    const fileName& dir,
    const fileType fType,
    const bool filter,
    const bool followLink
) const
{
    if (Pstream::parRun())
    {
        return uncollatedFileOperation::readDir
        (
            dir,
            fType,
            filter,
            followLink
        );
    }

    const fileName procDir(filePath(equivalentLagrangian(dir)));
    if (procDir.size())
    {
        const fileNameList contents
        (
            uncollatedFileOperation::readDir
            (
                procDir,
                fType,
                filter,
                followLink
            )
        );

        if (debug)
        {
            Pout<< indent
                << "autoReconstructFileOperation::readDir :"
                << " Returning from lagrangian directory searching:"
                << endl << indent
                << "    dir      :" << dir << endl << indent
                << "    procDir  :" << procDir << endl << indent
                << "    contents :";
            write(Pout, contents);
            Pout<< contents << endl << endl;
        }
        return contents;
    }
    else
    {
        return uncollatedFileOperation::readDir
        (
            dir,
            fType,
            filter,
            followLink
        );
    }

//    fileNameList contents;
//    if (dir.name() == cloud::prefix)
//    {
//        // Get equivalent processor directory
//
//        fileName pathAndInstance(dir.path());
//        fileName instance(pathAndInstance.name());
//        fileName procDir
//        (
//            filePath
//            (
//                pathAndInstance.path()/"processor0"/instance/cloud::prefix
//            )
//        );
//        if (procDir.size())
//        {
//            contents = uncollatedFileOperation::readDir
//            (
//                procDir,
//                fType,
//                filter,
//                followLink
//            );
//
//            if (debug)
//            {
//                Pout<< indent
//                    << "autoReconstructFileOperation::readDir :"
//                    << " Returning from lagrangian directory searching:"
//                    << endl << indent
//                    << "    procDir  :" << procDir
//                    << endl << indent
//                    << "    contents :" << contents << endl << endl;
//            }
//            return contents;
//        }
//    }
//
//    contents = uncollatedFileOperation::readDir
//    (
//        dir,
//        fType,
//        filter,
//        followLink
//    );
//
//    if (debug)
//    {
//        Pout<< indent
//            << "autoReconstructFileOperation::readDir :"
//            << " Returning from directory searching:"
//            << endl << indent
//            << "    dir      :" << dir
//            << endl << indent
//            << "    contents :" << contents << endl << endl;
//    }


//
//
////     // Force caching of processor directories. This makes sure
////     // that follow-on
////     // readDir does not need to do directory reading again
////     (void)lookupProcessorsPath(dir/"processor0");
//
//    Pout<< indent
//        << "autoReconstructFileOperation::readDir :"
//        << " dir:" << dir
//        << " fileType:" << label(fType)
//        << endl;
//
//    fileNameList contents;
//
//    // See if we have dir/processor0 (or dir/processorsXXX)
//    fileName procDir(filePath(dir/"processor0"));
//    if (procDir.size())
//    {
//        contents = uncollatedFileOperation::readDir
//        (
//            procDir,
//            fType,
//            filter,
//            followLink
//        );
//
//        if (debug)
//        {
//            Pout<< indent
//                << "autoReconstructFileOperation::readDir :"
//                << " Returning from (processor) directory searching:"
//                << endl << indent
//                << "    procDir  :" << procDir
//                << endl << indent
//                << "    contents :" << contents << endl << endl;
//        }
//    }
//    else
//    {
//        contents = uncollatedFileOperation::readDir
//        (
//            procDir,
//            fType,
//            filter,
//            followLink
//        );
//        if (debug)
//        {
//            Pout<< indent
//                << "autoReconstructFileOperation::readDir :"
//                << " Returning from directory searching:"
//                << endl << indent
//                << "    dir      :" << dir
//                << endl << indent
//                << "    contents :" << contents << endl << endl;
//        }
//    }
//
//
////     fileNameList contents
//     (
//         uncollatedFileOperation::readDir
//         (
//             dir,
//             fileType::directory,
//             filter,
//             followLink
//         )
//     );
//     Pout<< indent
//         << "autoReconstructFileOperation::readDir :"
//         << " contents:" << contents << endl;
//
//     fileName procDir;
//     forAll(contents, i)
//     {
//         const fileName& dirN = contents[i];
//
//         std::string::size_type pos = dirN.find("processor");
//         if (pos == 0)
//         {
//             procDir = dirN;
//             break;
//         }
//         else if (pos > 0 && dirN[pos-1] == '/')
//         {
//             procDir = dirN;
//             break;
//         }
//     }
//
// DebugVar(procDir);
//
//     if (procDir.empty())
//     {
//         if (fType != fileType::directory)
//         {
//             // Re-do
//             contents = uncollatedFileOperation::readDir
//             (
//                 dir,
//                 fType,
//                 filter,
//                 followLink
//             );
//         }
//     }
//     else
//     {
// DebugVar(dir/procDir);
//         contents = uncollatedFileOperation::readDir
//         (
//             dir/procDir,
//             fType,
//             filter,
//             followLink
//         );
//     }
//
//     if (debug)
//     {
//         Pout<< indent
//             << "autoReconstructFileOperation::readDir :"
//             << " Returning from directory searching:" << endl << indent
//             << "    dir      :" << dir << " fType:" << label(fType)
//             << endl << indent
//             << "    contents :" << contents << endl << endl;
//     }
//    return contents;
}


Foam::fileName Foam::fileOperations::autoReconstructFileOperation::dirPath
(
    const bool checkGlobal,
    const IOobject& io
) const
{
    fileName objPath;
    fileName procPath;
    if (haveProcPath(io, procPath))
    {
        objPath = io.objectPath();
    }
    else
    {
        // Try uncollated searching
        objPath = uncollatedFileOperation::filePath
        (
            checkGlobal,
            io,
            typeName
        );
    }

    if (debug)
    {
        Pout<< indent
            << "autoReconstructFileOperation::dirPath :"
            << " Returning from directory searching:" << endl << indent
            << "    objectPath:" << io.objectPath() << endl << indent
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileNameList
Foam::fileOperations::autoReconstructFileOperation::readObjects
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
        // tbd: lagrangian. dbDir="", local = "lagrangian/KinematicCloud"
        fileName path
        (
            filePath(db.path("processor0"/instance, db.dbDir()/local))
        );


        if (debug)
        {
            Pout<< indent
                << "autoReconstructFileOperation::readObjects :"
                << endl << indent
                << "    db       :" << db.path() << endl << indent
                << "    instance :" << instance << endl << indent
                << "    local    :" << local << endl << indent
                << "    store    :" << db.path(instance, db.dbDir()/local)
                << endl << indent
                << "    path     :" << path << endl;
        }


        if (path.size())    //Foam::isDir(path))
        {
            newInstance = instance;
            objects = Foam::readDir(path, fileType::file);

            if (debug)
            {
                Pout<< indent
                    << "autoReconstructFileOperation::readObjects :"
                    << " Returning processor directory searching:"
                    << endl << indent
                    << "    store:" << db.path(instance, db.dbDir()/local)
                    << endl << indent
                    << "    objects  :";
                write(Pout, objects);
                Pout<< endl << endl;
            }

            procObjects_.insert
            (
                db.path(instance, db.dbDir()/local),
                new fileNameList(objects)
            );

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

    if (debug)
    {
        Pout<< indent
            << "autoReconstructFileOperation::readObjects :"
            << " Returning from directory searching:" << endl << indent
            << "    path     :" << db.path(instance, db.dbDir()/local)
            << endl << indent
            << "    objects  :";
        write(Pout, objects);
        Pout<< endl << indent
            << "    newInst  :" << newInstance << endl << endl;
    }
    return objects;
}


Foam::instantList
Foam::fileOperations::autoReconstructFileOperation::findTimes
(
    const fileName& dir,
    const word& constantName
) const
{
    // Q: do we also include non-parallel directories?

    instantList times;//(uncollatedFileOperation::findTimes(dir, constantName));

    if (!Pstream::parRun())
    {
        fileName procDir(filePath(dir/"processor0"));
        if (procDir.size()) //exists(procDir))
        {
            instantList procTimes = uncollatedFileOperation::findTimes
            (
                procDir,
                constantName
            );
            times.append(procTimes);

            // Replicate uniqueOrder but with special handling of "constant"
            // (is not same as '0')
            const lessOrConstant compareOp(times);
            labelList order;
            sortedOrder(times, order, compareOp);
            if (order.size() > 1)
            {
                label n = 0;
                for (label i = 0; i < order.size() - 1; ++i)
                {
                    if (!compareOp.equal(order[i], order[i+1]))
                    {
                        order[n++] = order[i];
                    }
                }
                order[n++] = order[order.size()-1];
                order.setSize(n);
            }
            times = UIndirectList<instant>(times, order)();
            if (debug)
            {
                Pout<< indent
                    << "autoReconstructFileOperation::findTimes :"
                    << " Returning from processor time searching:" << endl
                    << indent
                    << "    procDir :" << dir << endl << indent
                    << "    times   :";
                write(Pout, times);
                Pout<< endl << endl;
            }
        }
        else
        {
            times = uncollatedFileOperation::findTimes(dir, constantName);
            if (debug)
            {
                Pout<< indent
                    << "autoReconstructFileOperation::findTimes :"
                    << " Returning from time searching:" << endl << indent
                    << "    dir   :" << dir << endl << indent
                    << "    times :";
                write(Pout, times);
                Pout<< endl << endl;
            }
        }
    }
    else
    {
        times = uncollatedFileOperation::findTimes(dir, constantName);
        if (debug)
        {
            Pout<< indent
                << "autoReconstructFileOperation::findTimes :"
                << " Returning from time searching:" << endl << indent
                << "    dir   :" << dir << endl << indent
                << "    times :";
            write(Pout, times);
            Pout<< endl << endl;
        }
    }
    return times;
}


// bool Foam::fileOperations::autoReconstructFileOperation::readHeader
// (
//     IOobject& io,
//     const fileName& fName,
//     const word& typeName
// ) const
// {
//     if (debug)
//     {
//         Pout<< indent
//             << "autoReconstructFileOperation::readHeader :"
//             << " Reading:" << typeName << " from: " << fName << endl;
//     }
//     bool ok = uncollatedFileOperation::readHeader(io, fName, typeName);
//
//     if (debug)
//     {
//         Pout<< indent
//             << "autoReconstructFileOperation::readHeader :" << endl
//             << indent
//             << "    fName :" << fName << endl << indent
//             << "    ok    :" << ok << endl << endl;
//     }
//     return ok;
// }
//
//
//
Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::autoReconstructFileOperation::readStream
(
    regIOobject& io,
    const fileName& fName,
    const word& type,
    const bool valid
) const
{
    if (debug)
    {
        Pout<< indent
            << "autoReconstructFileOperation::readStream :"
            << endl << indent
            << "    io    :" << io.objectPath() << endl << indent
            << "    local    :" << io.local() << endl << indent
            << "    fName :" << fName << endl << indent
            << "    type  :" << type << endl << endl;
    }

    autoPtr<ISstream> isPtr;

    if (!valid)
    {
        isPtr = autoPtr<ISstream>(new dummyISstream());
        return isPtr;
    }

    autoPtr<streamReconstructor> reconstructor
    (
        streamReconstructor::New(type)
    );

    if (reconstructor.valid())
    {
        if (debug)
        {
            Pout<< indent
                << "autoReconstructFileOperation::readStream :"
                << " Found reconstructor for type:" << type
                << " of object: " << io.objectPath() << endl;
        }

        OStringStream os(IOstream::BINARY);

        // Save current file handler and install basic ('uncollated') one
Pout<< "**INSTALLING BASIC**:" << basicFileHandler_().type() << endl;
        autoPtr<fileOperation> orig(fileOperation::fileHandlerPtr_.ptr());
        (void)fileHandler(basicFileHandler_);

        bool ok = reconstructor().reconstruct(io, false, os);

        // Restore current file handler
Pout<< "**UNINSTALLING BASIC**:" << fileHandler().type() << endl;
        basicFileHandler_ = fileOperation::fileHandlerPtr_.ptr();
        (void)fileHandler(orig);

        if (!ok)
        {
            //isPtr.reset(dummyISstream());
        }
        else
        {
            isPtr.reset(new IStringStream(os.str(), IOstream::BINARY));
        }
    }
    else
    {
        isPtr = uncollatedFileOperation::readStream(io, fName, type, valid);
    }

    if (debug)
    {
        Pout<< indent
            << "autoReconstructFileOperation::readStream :"
            << endl << indent
            << "    fName :" << fName << endl << indent
            << "    isPtr :" << isPtr.valid() << endl << endl;
    }
    return isPtr;
}


bool Foam::fileOperations::autoReconstructFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& type
) const
{
    if (debug)
    {
        Pout<< indent
            << "autoReconstructFileOperation::read :"
            << " Reading:" << type << " from: " << io.objectPath() << endl;

        Pout<< indent << "    local    :" << io.local() << endl;
    }

    fileName procPath;
    if (haveProcPath(io, procPath))
    {
        if (debug)
        {
            Pout<< indent
                << "autoReconstructFileOperation::read :"
                << " Searching for reconstructor for type:" << type
                << " of object: " << io.objectPath() << endl;
        }
        autoPtr<streamReconstructor> reconstructor
        (
            streamReconstructor::New(type)
        );

        if (reconstructor.valid())
        {
            if (debug)
            {
                Pout<< indent
                    << "autoReconstructFileOperation::read :"
                    << " Found reconstructor for type:" << type
                    << " of object: " << io.objectPath() << endl;
            }

            OStringStream os(IOstream::BINARY);

            // Save current file handler and install basic ('uncollated') one
Pout<< "**INSTALLING BASIC**:" << basicFileHandler_().type() << endl;
            autoPtr<fileOperation> orig(fileOperation::fileHandlerPtr_.ptr());
            (void)fileHandler(basicFileHandler_);

            bool ok = reconstructor().reconstruct(io, false, os);

            // Restore current file handler
Pout<< "**UNINSTALLING BASIC**:" << fileHandler().type() << endl;
            basicFileHandler_ = fileOperation::fileHandlerPtr_.ptr();
            (void)fileHandler(orig);

            if (!ok)
            {
                return false;
            }
            else
            {
                IStringStream is(os.str(), IOstream::BINARY);
                return io.readData(is);
            }
        }
        else
        {
            return uncollatedFileOperation::read(io, masterOnly, format, type);
        }
    }
    else
    {
        return uncollatedFileOperation::read(io, masterOnly, format, type);
    }
}


// ************************************************************************* //
