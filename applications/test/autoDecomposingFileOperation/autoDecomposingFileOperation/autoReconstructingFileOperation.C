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

#include "autoReconstructingFileOperation.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvFieldDecomposer.H"
#include "addToRunTimeSelectionTable.H"
#include "processorMeshes.H"
#include "fvFieldReconstructor.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(autoReconstructingFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        autoReconstructingFileOperation,
        word
    );


//     class installFileOp
//     {
//     public:
//
//         installFileOp()
//         {
//             // Install autoReconstructing as fileHandler
//             autoPtr<fileOperation> handler
//             (
//                 new autoReconstructingFileOperation(true)
//             );
//             Foam::fileHandler(handler);
//         }
//
//         ~installFileOp()
//         {
//             autoPtr<fileOperation> handler(nullptr);
//             Foam::fileHandler(handler);
//         }
//     };
//     installFileOp installFileOp_;
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileOperations::autoReconstructingFileOperation::haveProcPath
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
        procObjectPath = fileName
        (
            io.rootPath()
           /io.caseName()
           /"processor0"
           /io.instance()
           /io.db().dbDir()
           /io.local()
           /io.name()
        );
        return exists(procObjectPath);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::autoReconstructingFileOperation::
autoReconstructingFileOperation
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

Foam::fileOperations::autoReconstructingFileOperation::
~autoReconstructingFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::autoReconstructingFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& typeName
) const
{
    if (debug)
    {
        Pout<< "autoReconstructingFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

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
        Pout<< "autoReconstructingFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileName Foam::fileOperations::autoReconstructingFileOperation::dirPath
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
        Pout<< "autoReconstructingFileOperation::dirPath :"
            << " Returning from directory searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileNameList
Foam::fileOperations::autoReconstructingFileOperation::readObjects
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
        fileName path(db.path("processor0"/instance, db.dbDir()/local));

DebugVar(path);

        if (Foam::isDir(path))
        {
            newInstance = instance;
            objects = Foam::readDir(path, fileName::FILE);
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
        Pout<< "autoReconstructingFileOperation::readObjects :"
            << " Returning from directory searching:" << endl
            << "    path     :" << db.path() << endl
            << "    objects  :" << objects << endl << endl;
    }
    return objects;
}


Foam::instantList
Foam::fileOperations::autoReconstructingFileOperation::findTimes
(
    const fileName& dir,
    const word& constantName
) const
{
    instantList times(uncollatedFileOperation::findTimes(dir, constantName));

    if (!Pstream::parRun())
    {
        fileName procDir(dir/"processor0");
        if (exists(procDir))
        {
            instantList procTimes = uncollatedFileOperation::findTimes
            (
                procDir,
                constantName
            );
            bool procHasConstant =
            (
                procTimes.size() > 0
             && procTimes[0].name() == constantName
            );

            if (procHasConstant)
            {
                SubList<instant> realTimes(procTimes, procTimes.size()-1, 1);
                times.append(realTimes);
            }
            else
            {
                times.append(procTimes);
            }

            bool hasConstant =
            (
                times.size() > 0
             && times[0].name() == constantName
            );
            if (hasConstant)
            {
                SubList<instant> realTimes(times, times.size()-1, 1);
                labelList order;
                uniqueOrder(realTimes, order);

                instantList newTimes(order.size()+1);
                newTimes[0] = times[0];
                forAll(order, i)
                {
                    newTimes[i+1] = realTimes[order[i]];
                }
                times.transfer(newTimes);
            }
            else
            {
                labelList order;
                uniqueOrder(times, order);
                times = UIndirectList<instant>(times, order)();
            }
        }
    }

    if (debug)
    {
        Pout<< "autoReconstructingFileOperation::findTimes :"
            << " Returning from time searching:" << endl
            << "    dir   :" << dir << endl
            << "    times :" << times << endl << endl;
    }
    return times;
}


bool Foam::fileOperations::autoReconstructingFileOperation::readHeader
(
    IOobject& io,
    const fileName& fName,
    const word& typeName
) const
{
    //if (debug)
    //{
    //    Pout<< "autoReconstructingFileOperation::readHeader :"
    //        << " Reading:" << type << " from: " << fName << endl;
    //}
    if (typeName == volScalarField::typeName)
    {
    }

    bool ok = uncollatedFileOperation::readHeader(io, fName, typeName);

    if (debug)
    {
        Pout<< "autoReconstructingFileOperation::readHeader :" << endl
            << "    fName :" << fName << endl
            << "    ok    :" << ok << endl << endl;
    }
    return ok;
}



// Foam::autoPtr<Foam::ISstream>
// Foam::fileOperations::autoReconstructingFileOperation::readStream
// (
//     regIOobject& io,
//     const fileName& fName,
//     const word& typeName,
//     const bool valid
// ) const
// {
//     autoPtr<ISstream> isPtr
//     (
//         uncollatedFileOperation::readStream
//         (
//             io,
//             fName,
//             typeName,
//             valid
//         )
//     );
//
//     if (debug)
//     {
//         Pout<< "autoReconstructingFileOperation::readStream :" << endl
//             << "    fName :" << fName << endl
//             << "    isPtr :" << isPtr.valid() << endl << endl;
//     }
//     return isPtr;
// }


bool Foam::fileOperations::autoReconstructingFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& type
) const
{
    if (debug)
    {
        Pout<< "autoReconstructingFileOperation::read :"
            << " Reading:" << type << " from: " << io.objectPath() << endl;
    }

    fileName procPath;
    if (haveProcPath(io, procPath))
    {
        DebugVar(procPath);
    }
    return uncollatedFileOperation::read
    (
        io,
        masterOnly,
        format,
        type
    );
}


// ************************************************************************* //
