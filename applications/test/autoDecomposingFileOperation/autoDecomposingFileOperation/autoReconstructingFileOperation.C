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
        Pout<< indent
            << "autoReconstructingFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    fileName objPath;
    fileName procPath;
    if (haveProcPath(io, procPath))
    {
        objPath = procPath; //io.objectPath();

        if (debug)
        {
            Pout<< indent
                << "autoReconstructingFileOperation::filePath :"
                << " Returning from processor searching:" << endl << indent
                << "    objectPath:" << io.objectPath() << " type:" << typeName
                << endl << indent
                << "    filePath  :" << objPath << endl << endl;
        }
        return objPath;
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
            << "autoReconstructingFileOperation::filePath :"
            << " Returning from file searching:" << endl << indent
            << "    objectPath:" << io.objectPath() << " type:" << typeName
            << endl << indent
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
        Pout<< indent
            << "autoReconstructingFileOperation::dirPath :"
            << " Returning from directory searching:" << endl << indent
            << "    objectPath:" << io.objectPath() << endl << indent
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

            if (debug)
            {
                Pout<< indent
                    << "autoReconstructingFileOperation::readObjects :"
                    << " Returning processor directory searching:"
                    << endl << indent
                    << "    path     :" << db.path() << endl << indent
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
            << "autoReconstructingFileOperation::readObjects :"
            << " Returning from directory searching:" << endl << indent
            << "    path     :" << db.path() << endl << indent
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

//     if (debug)
//     {
//         Pout<< indent
//             << "autoReconstructingFileOperation::findTimes :"
//             << " Returning from time searching:" << endl << indent
//             << "    dir   :" << dir << endl << indent
//             << "    times :" << times << endl << endl;
//     }
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
    //    Pout<< indent
    //        << "autoReconstructingFileOperation::readHeader :"
    //        << " Reading:" << type << " from: " << fName << endl;
    //}
    if (typeName == volScalarField::typeName)
    {
    }

    bool ok = uncollatedFileOperation::readHeader(io, fName, typeName);

    if (debug)
    {
        Pout<< indent
            << "autoReconstructingFileOperation::readHeader :" << endl << indent
            << "    fName :" << fName << endl << indent
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
//         Pout<< indent
//             << "autoReconstructingFileOperation::readStream :"
//             << endl << indent
//             << "    fName :" << fName << endl << indent
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
        Pout<< indent
            << "autoReconstructingFileOperation::read :"
            << " Reading:" << type << " from: " << io.objectPath() << endl;
    }

    fileName procPath;
    if (haveProcPath(io, procPath) && type == volScalarField::typeName)
    {
        DebugVar(procPath);
        // Read file from procPath and reconstruct instead of from objectPath

        const label nProcs = 2;

        // Create the processor databases
        PtrList<Time> databases(nProcs);

        forAll(databases, proci)
        {
            Pout<< indent << "** Loading Time for:" << proci << nl << endl;
            Pout<< incrIndent;

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
            databases[proci].setTime(io.time());

            Pout<< decrIndent;
            Pout<< indent << "** DONE Loading Time for:" << proci << nl << endl;
        }

        Pout<< indent << "** Loading processors" << nl << endl;
        Pout<< incrIndent;
        // Read all meshes and addressing to reconstructed mesh
        processorMeshes procMeshes(databases, polyMesh::defaultRegion);
        Pout<< decrIndent;
        Pout<< indent << "** DONE Loading processors" << nl << endl;

DebugVar(io.objectPath());
DebugVar(io.db().name());
DebugVar(io.db().path());

DebugVar(procMeshes.boundaryProcAddressing());
DebugVar(procMeshes.cellProcAddressing());
DebugVar(procMeshes.faceProcAddressing());

        fvMesh& mesh = const_cast<fvMesh&>
        (
            dynamic_cast<const fvMesh&>(io.db())
        );
        Pout<< "mesh nCells:" << mesh.nCells() << endl;
        Pout<< "mesh path:" << mesh.polyMesh::path() << endl;

        // Construct reconstructor
        fvFieldReconstructor fvReconstructor
        (
            const_cast<fvMesh&>(dynamic_cast<const fvMesh&>(io.db())),
            procMeshes.meshes(),
            procMeshes.faceProcAddressing(),
            procMeshes.cellProcAddressing(),
            procMeshes.boundaryProcAddressing()
        );

        // Read all fields
        PtrList<volScalarField> procFields(nProcs);
        {
            for (label proci = 0; proci < nProcs; proci++)
            {
                const fvMesh& procMesh = procMeshes.meshes()[proci];

                IOobject procIO(io, procMesh);

Pout<< indent << "For proc:" << proci << " trying to read "
    << procIO.objectPath() << endl;
                Pout<< incrIndent;
                procFields.set(proci, new volScalarField(procIO, procMesh));
                Pout<< decrIndent;
            }
        }

        volScalarField& ioFld = dynamic_cast<volScalarField&>(io);

DebugVar(procFields);
        tmp<volScalarField> tfld
        (
            fvReconstructor.reconstructFvVolumeField
            (
                io,
                procFields
            )
        );
        ioFld = tfld;
        return true;
    }
    return uncollatedFileOperation::read(io, masterOnly, format, type);
}


// ************************************************************************* //
