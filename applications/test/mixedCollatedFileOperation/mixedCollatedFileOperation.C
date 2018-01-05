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

#include "mixedCollatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "threadedCollatedOFstream.H"
#include "decomposedBlockData.H"
#include "registerSwitch.H"
#include "masterOFstream.H"
#include "OFstream.H"
#include "collatedFileOperation.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(mixedCollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        mixedCollatedFileOperation,
        word
    );

    float mixedCollatedFileOperation::maxThreadFileBufferSize
    (
        debug::floatOptimisationSwitch("maxThreadFileBufferSize", 1e9)
    );
    registerOptSwitch
    (
        "maxThreadFileBufferSize",
        float,
        mixedCollatedFileOperation::maxThreadFileBufferSize
    );

    // Mark as needing threaded mpi
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        collatedFileOperationInitialise,
        word,
        mixedCollated
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// bool Foam::fileOperations::mixedCollatedFileOperation::appendObject
// (
//     const regIOobject& io,
//     const fileName& pathName,
//     IOstream::streamFormat fmt,
//     IOstream::versionNumber ver,
//     IOstream::compressionType cmp
// ) const
// {
//     // Append to processors/ file
//
//     label proci = detectProcessorPath(io.objectPath());
//
//     if (debug)
//     {
//         Pout<< "writeObject:" << " : For local object : "
//             << io.name()
//             << " appending processor " << proci
//             << " isMaster:" << isMaster(nProcs_, proci)
//             << " data to " << pathName << endl;
//     }
//
//     if (proci == -1)
//     {
//         FatalErrorInFunction
//             << "Not a valid processor path " << pathName
//             << exit(FatalError);
//     }
//
//     bool master = isMaster(nProcs_, proci);
//
//
//     // Create string from all data to write
//     string buf;
//     {
//         OStringStream os(fmt, ver);
//         if (master)
//         {
//             if (!io.writeHeader(os))
//             {
//                 return false;
//             }
//         }
//
//         // Write the data to the Ostream
//         if (!io.writeData(os))
//         {
//             return false;
//         }
//
//         if (master)
//         {
//             IOobject::writeEndDivider(os);
//         }
//
//         buf = os.str();
//     }
//
//
//     // Note: cannot do append + compression. This is a limitation
//     // of ogzstream (or rather most compressed formats)
//
//     OFstream os
//     (
//         pathName,
//         IOstream::BINARY,
//         ver,
//         IOstream::UNCOMPRESSED, // no compression
//         !master
//     );
//
//     if (!os.good())
//     {
//         FatalIOErrorInFunction(os)
//             << "Cannot open for appending"
//             << exit(FatalIOError);
//     }
//
//     if (master)
//     {
//         IOobject::writeBanner(os)
//             << "FoamFile\n{\n"
//             << "    version     " << os.version() << ";\n"
//             << "    format      " << os.format() << ";\n"
//             << "    class       " << decomposedBlockData::typeName
//             << ";\n"
//             << "    location    " << pathName << ";\n"
//             << "    object      " << pathName.name() << ";\n"
//             << "}" << nl;
//         IOobject::writeDivider(os) << nl;
//     }
//
//     // Write data
//     UList<char> slice
//     (
//         const_cast<char*>(buf.data()),
//         label(buf.size())
//     );
//     os << nl << "// Processor" << proci << nl << slice << nl;
//
//     return os.good();
// }


Foam::labelPair Foam::fileOperations::mixedCollatedFileOperation::commsGroup
(
    const label nProcs,
    const label proci
)
{
    label n = 2;    //1;
    label master = (proci / n) * n;
    return labelPair(master, n);
}


bool Foam::fileOperations::mixedCollatedFileOperation::isMaster
(
    const label nProcs,
    const label proci
)
{
    return commsGroup(nProcs, proci).first() == proci;
}


Foam::labelList Foam::fileOperations::mixedCollatedFileOperation::subRanks
(
    const label n
)
{
//     const string myHostName(hostName());
//
//     stringList hosts(Pstream::nProcs());
//     hosts[Pstream::myProcNo()] = myHostName;
//     Pstream::gatherList(hosts);
//     Pstream::scatterList(hosts);
//
//     // Collect procs with same hostname
//     DynamicList<label> subRanks(64);
//     forAll(hosts, proci)
//     {
//         if (hosts[proci] == myHostName)
//         {
//             subRanks.append(proci);
//         }
//     }

    labelPair group(commsGroup(n, Pstream::myProcNo()));
    label masterProc = group.first();
    label nProcs = group.second();

    labelList subRanks(nProcs);
    forAll(subRanks, i)
    {
        subRanks[i] = masterProc++;
    }

DebugVar(subRanks);
    return subRanks;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::mixedCollatedFileOperation::mixedCollatedFileOperation
(
    const bool verbose
)
:
    //masterUncollatedFileOperation
    collatedFileOperation
    (
        UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            subRanks(Pstream::nProcs())
        ),
        false
    )
//    writer_(maxThreadFileBufferSize, comm_),
//    nProcs_(-1)
{
//     const List<int>& procs(UPstream::procID(comm_));
//
//     if (Pstream::parRun())
//     {
//         processorsDir_ =
//             processorsBaseDir
//           + Foam::name(Pstream::nProcs())
//           + "_"
//           + Foam::name(procs[0])
//           + "-"
//           + Foam::name(procs.last());
//     }
//     else
//     {
//         processorsDir_ = processorsBaseDir;
//     }
//
//     if (verbose)
//     {
//         Info<< "I/O    : " << typeName
//             << " (output directory " << processorsDir_
//             << ", maxThreadFileBufferSize " << maxThreadFileBufferSize
//             << ')' << endl;
//
//         if (maxThreadFileBufferSize == 0)
//         {
//             Info<< "         Threading not activated "
//                    "since maxThreadFileBufferSize = 0." << nl
//                 << "         Writing may run slowly for large file sizes."
//                 << endl;
//         }
//         else
//         {
//             Info<< "         Threading activated "
//                    "since maxThreadFileBufferSize > 0." << nl
//                 << "         Requires large enough buffer to collect all data"
//                     " or thread support " << nl
//                 << "         enabled in MPI. If thread support cannot be "
//                    "enabled, deactivate" << nl
//                 << "         threading by setting maxThreadFileBufferSize "
//                     "to 0 in" << nl
//                 << "         $FOAM_ETC/controlDict"
//                 << endl;
//         }
//
//         if
//         (
//             regIOobject::fileModificationChecking
//          == regIOobject::inotifyMaster
//         )
//         {
//             WarningInFunction
//                 << "Resetting fileModificationChecking to inotify" << endl;
//         }
//
//         if
//         (
//             regIOobject::fileModificationChecking
//          == regIOobject::timeStampMaster
//         )
//         {
//             WarningInFunction
//                 << "Resetting fileModificationChecking to timeStamp" << endl;
//         }
//     }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::mixedCollatedFileOperation::~mixedCollatedFileOperation()
{
    if (comm_ != -1)
    {
        UPstream::freeCommunicator(comm_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Foam::fileName Foam::fileOperations::mixedCollatedFileOperation::objectPath
// (
//     const IOobject& io,
//     const word& typeName
// ) const
// {
//     // Replacement for objectPath
//     if (io.time().processorCase())
//     {
//         return masterUncollatedFileOperation::objectPath
//         (
//             io,
//             fileOperation::PROCESSORSOBJECT,
//             io.instance()
//         );
//     }
//     else
//     {
//         return masterUncollatedFileOperation::objectPath
//         (
//             io,
//             fileOperation::OBJECT,
//             io.instance()
//         );
//     }
// }
//
//
// bool Foam::fileOperations::mixedCollatedFileOperation::writeObject
// (
//     const regIOobject& io,
//     IOstream::streamFormat fmt,
//     IOstream::versionNumber ver,
//     IOstream::compressionType cmp,
//     const bool valid
// ) const
// {
//     const Time& tm = io.time();
//     const fileName& inst = io.instance();
//
//     if (inst.isAbsolute() || !tm.processorCase())
//     {
// DebugVar(io.path());
//         mkDir(io.path());
//         fileName pathName(io.objectPath());
//
//         if (debug)
//         {
//             Pout<< "writeObject:"
//                 << " : For object : " << io.name()
//                 << " falling back to master-only output to " << io.path()
//                 << endl;
//         }
//
//         masterOFstream os
//         (
//             pathName,
//             fmt,
//             ver,
//             cmp,
//             false,
//             valid
//         );
//
//         // If any of these fail, return (leave error handling to Ostream class)
//         if (!os.good())
//         {
//             return false;
//         }
//         if (!io.writeHeader(os))
//         {
//             return false;
//         }
//         // Write the data to the Ostream
//         if (!io.writeData(os))
//         {
//             return false;
//         }
//         IOobject::writeEndDivider(os);
//
//         return true;
//     }
//     else
//     {
//         if (io.global())
//         {
//             //const word actualProcsDir(processorsDir(io));
//
//             // Construct the equivalent processors/ directory
//             fileName path(processorsPath(io, inst, processorsDir(io)));
//
// DebugVar(path);
//             mkDir(path);
//             fileName pathName(path/io.name());
//
//             if (debug)
//             {
//                 Pout<< "writeObject:" << " : For global object : " << io.name()
//                     << " falling back to master-only output to " << pathName
//                     << endl;
//             }
//
//             masterOFstream os
//             (
//                 pathName,
//                 fmt,
//                 ver,
//                 cmp,
//                 false,
//                 valid
//             );
//
//             // If any of these fail, return (leave error handling to Ostream
//             // class)
//             if (!os.good())
//             {
//                 return false;
//             }
//             if (!io.writeHeader(os))
//             {
//                 return false;
//             }
//             // Write the data to the Ostream
//             if (!io.writeData(os))
//             {
//                 return false;
//             }
//             IOobject::writeEndDivider(os);
//
//             return true;
//         }
//         else if (!Pstream::parRun())
//         {
//             // Special path for e.g. decomposePar. Append to
//             // processorsDDD/ file
//
//             // Re-construct the equivalent processors/ directory
//             fileName path(processorsPath(io, inst, processorsDir(io)));
//             mkDir(path);
//             fileName pathName(path/io.name());
//
//             if (debug)
//             {
//                 Pout<< "writeObject:"
//                     << " : For object : " << io.name()
//                     << " appending to " << pathName << endl;
//             }
//
//             return appendObject(io, pathName, fmt, ver, cmp);
//         }
//         else
//         {
//             // Construct the equivalent processors/ directory
//             fileName path(processorsPath(io, inst, processorsDir(io)));
//
//             mkDir(path);
//             fileName pathName(path/io.name());
//
//             if (debug)
//             {
//                 Pout<< "writeObject:"
//                     << " : For object : " << io.name()
//                     << " starting collating output to " << pathName << endl;
//             }
//
//             threadedCollatedOFstream os(writer_, pathName, fmt, ver, cmp);
//
//             // If any of these fail, return (leave error handling to Ostream
//             // class)
//             if (!os.good())
//             {
//                 return false;
//             }
//             if (Pstream::master(comm_) && !io.writeHeader(os))
//             {
//                 return false;
//             }
//             // Write the data to the Ostream
//             if (!io.writeData(os))
//             {
//                 return false;
//             }
//             if (Pstream::master(comm_))
//             {
//                 IOobject::writeEndDivider(os);
//             }
//
//             return true;
//         }
//     }
// }
//
//
// Foam::word Foam::fileOperations::mixedCollatedFileOperation::processorsDir
// (
//     const IOobject& io
// ) const
// {
//     if (Pstream::parRun())
//     {
//         // Use pre-calculated directory
//         return processorsDir_;
//     }
//     else
//     {
//         return processorsDir(io.objectPath());
//     }
// }
//
//
// Foam::word Foam::fileOperations::mixedCollatedFileOperation::processorsDir
// (
//     const fileName& fName
// ) const
// {
//     if (Pstream::parRun())
//     {
//         // Use pre-calculated directory
//         return processorsDir_;
//     }
//     else
//     {
//         // Detect processor number from the filename
//         label proci = detectProcessorPath(fName);
//
//         labelPair group(commsGroup(nProcs_, proci));
//         label groupMaster = group.first();
//         label groupSize = group.second();
//
//         word processorsDir
//         (
//             processorsBaseDir
//           + Foam::name(nProcs_)
//           + "_"
//           + Foam::name(groupMaster)
//           + "-"
//           + Foam::name(groupMaster+groupSize-1)
//         );
//
//         return processorsDir;
//     }
// }
//
//
// void Foam::fileOperations::mixedCollatedFileOperation::setNProcs
// (
//     const label nProcs
// )
// {
//     // Changed number of decompositions
//     nProcs_ = nProcs;
//
//     if (debug)
//     {
//         Pout<< "setNProcs:" << " : Setting nProcs to " << nProcs_ << endl;
//     }
// }


// ************************************************************************* //
