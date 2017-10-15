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

#include "OFstreamCollator.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "IOstreams.H"
#include "Pstream.H"
#include "decomposedBlockData.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstreamCollator, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::OFstreamCollator::writeFile
(
    const label comm,
    const word& typeName,
    const fileName& fName,
    const string& s,
    const List<char>& slaveData,
    const List<int>& slaveOffsets,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool append
)
{
    if (debug)
    {
        Pout<< "OFstreamCollator : Writing " << s.size()
            << " bytes to " << fName
            << " using comm " << comm << endl;
    }

    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm))
    {
        Foam::mkDir(fName.path());
        osPtr.reset
        (
            new OFstream
            (
                fName,
                fmt,
                ver,
                cmp,
                append
            )
        );

        // We don't have IOobject so cannot use IOobject::writeHeader
        OSstream& os = osPtr();
        decomposedBlockData::writeHeader
        (
            os,
            ver,
            fmt,
            typeName,
            "",
            fName,
            fName.name()
        );
    }

    UList<char> slice(const_cast<char*>(s.data()), label(s.size()));

    // Determine sizes to receive
    labelList recvSizes(Pstream::nProcs(comm));
    {
        char* data = reinterpret_cast<char*>(recvSizes.begin());

        List<int> recvOffsets(recvSizes.size());
        forAll(recvOffsets, proci)
        {
            recvOffsets[proci] =
                reinterpret_cast<char*>(&recvSizes[proci])
              - data;
        }
        label sz = slice.byteSize();
        UPstream::gather
        (
            reinterpret_cast<char*>(&sz),
            sizeof(sz),
            data,
            List<int>(recvSizes.size(), sizeof(label)),
            recvOffsets,
            comm
        );
    }

    // Assuming threaded writing hides any slowness so we
    // can use scheduled communication to send the data to
    // the master processor in order. However can be unstable so
    // default is non-blocking.

    List<std::streamoff> start;
    decomposedBlockData::writeBlocks
    (
        comm,
        osPtr,
        start,
        slice,
        slaveData,
        slaveOffsets,
        recvSizes,
        UPstream::commsTypes::nonBlocking,  //scheduled,
        false       // do not reduce return state
    );

    if (osPtr.valid() && !osPtr().good())
    {
        FatalIOErrorInFunction(osPtr())
            << "Failed writing to " << fName << exit(FatalIOError);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Finished writing " << s.size() << " bytes";
        if (UPstream::master(comm))
        {
            label sum = 0;
            forAll(recvSizes, i)
            {
                sum += recvSizes[i];
            }
            Pout<< " (overall " << sum << ")";
        }
        Pout<< " to " << fName
            << " using comm " << comm << endl;
    }

    return true;
}


void* Foam::OFstreamCollator::writeAll(void *threadarg)
{
    OFstreamCollator& handler = *static_cast<OFstreamCollator*>(threadarg);

    // Consume stack
    while (true)
    {
        writeData* ptr = nullptr;

        lockMutex(handler.mutex_);
        if (handler.objects_.size())
        {
            ptr = handler.objects_.pop();
        }
        unlockMutex(handler.mutex_);

        if (!ptr)
        {
            break;
        }
        else
        {
            bool ok = writeFile
            (
                handler.comm_,
                ptr->typeName_,
                ptr->pathName_,
                ptr->data_,
                ptr->format_,
                ptr->version_,
                ptr->compression_,
                ptr->append_
            );
            if (!ok)
            {
                FatalIOErrorInFunction(ptr->pathName_)
                    << "Failed writing " << ptr->pathName_
                    << exit(FatalIOError);
            }

            delete ptr;
        }
        //sleep(1);
    }

    if (debug)
    {
        Pout<< "OFstreamCollator : Exiting write thread " << endl;
    }

    lockMutex(handler.mutex_);
    handler.threadRunning_ = false;
    unlockMutex(handler.mutex_);

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstreamCollator::OFstreamCollator(const off_t maxBufferSize)
:
    maxBufferSize_(maxBufferSize),
    mutex_
    (
        maxBufferSize_ > 0
      ? allocateMutex()
      : -1
    ),
    thread_
    (
        maxBufferSize_ > 0
      ? allocateThread()
      : -1
    ),
    threadRunning_(false),
    comm_
    (
        UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            identity(UPstream::nProcs(UPstream::worldComm))
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstreamCollator::~OFstreamCollator()
{
    if (threadRunning_)
    {
        if (debug)
        {
            Pout<< "~OFstreamCollator : Waiting for write thread" << endl;
        }

        joinThread(thread_);
    }
    if (thread_ != -1)
    {
        freeThread(thread_);
    }
    if (mutex_ != -1)
    {
        freeMutex(mutex_);
    }
    if (comm_ != -1)
    {
        UPstream::freeCommunicator(comm_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OFstreamCollator::write
(
    const word& typeName,
    const fileName& fName,
    const string& data,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool append
)
{
    // Determine sizes to receive
    labelList recvSizes(Pstream::nProcs());
    {
        char* data = reinterpret_cast<char*>(recvSizes.begin());

        List<int> recvOffsets(recvSizes.size());
        forAll(recvOffsets, proci)
        {
            recvOffsets[proci] =
                reinterpret_cast<char*>(&recvSizes[proci])
              - data;
        }
        label sz = slice.byteSize();
        UPstream::gather
        (
            reinterpret_cast<char*>(&sz),
            sizeof(sz),
            data,
            List<int>(recvSizes.size(), sizeof(label)),
            recvOffsets,
            Pstream::worldComm      // Do NOT use thread communicator
        );
    }

DebugVar(recvSizes);
    label totalSlaveSize = sum
    (
        SubList<label>(recvSizes, recvSizes.size()-1, 1)
    );
DebugVar(totalSlaveSize);

    if
    (
        totalSlaveSize
      < fileOperations::masterUncollatedFileOperation::maxMasterFileBufferSize
    )
    {
        // Allocate local buffer for all collated data
        autoPtr<collatedData> fileAndDataPtr
        (
            new writeData(typeName, fName, data, fmt, ver, cmp, append)
        );
        collatedData& fileAndData = fileAndDataPtr();

        // Gather the slave data and insert into fileAndData
        UList<char> slice(const_cast<char*>(data.data()), label(data.size()));
        decomposedBlockData::gatherSlaveData
        (
            Pstream::worldComm,
            slice.begin(),
            recvSizes,

            1,                          // startProc,
            Pstream::nProcs()-1,        // n procs

            fileAndData.slaveOffsets_,
            fileAndData.data_
        );

DebugVar(fileAndData.slaveOffsets_);

        // Append to thread buffer
        lockMutex(mutex_);
        objects_.push(fileAndDataPtr.ptr());
        unlockMutex(mutex_);

        // Start thread if not running
        lockMutex(mutex_);
        if (!threadRunning_)
        {
            createThread(thread_, writeAll, this);
            if (debug)
            {
                Pout<< "OFstreamCollator : Started write thread " << endl;
            }
            threadRunning_ = true;
        }
        unlockMutex(mutex_);
    }
    else if (maxBufferSize_ > 0)
    {
        while (true)
        {
            // Count files to be written
            off_t totalSize = 0;

            lockMutex(mutex_);
            forAllConstIter(FIFOStack<writeData*>, objects_, iter)
            {
                totalSize += iter()->size();
            }
            unlockMutex(mutex_);

            if
            (
                totalSize == 0
             || (totalSize+off_t(data.size()) < maxBufferSize_)
            )
            {
                break;
            }

            if (debug)
            {
                lockMutex(mutex_);
                Pout<< "OFstreamCollator : Waiting for buffer space."
                    << " Currently in use:" << totalSize
                    << " limit:" << maxBufferSize_
                    << " files:" << objects_.size()
                    << endl;
                unlockMutex(mutex_);
            }

            sleep(5);
        }

        if (debug)
        {
            Pout<< "OFstreamCollator : relaying write of " << fName
                << " to thread " << endl;
        }

        lockMutex(mutex_);
        objects_.push
        (
            new writeData(typeName, fName, data, fmt, ver, cmp, append)
        );
        unlockMutex(mutex_);

        lockMutex(mutex_);
        if (!threadRunning_)
        {
            createThread(thread_, writeAll, this);
            if (debug)
            {
                Pout<< "OFstreamCollator : Started write thread " << endl;
            }
            threadRunning_ = true;
        }
        unlockMutex(mutex_);

        return true;
    }
    else
    {
        // Immediate writing
        return writeFile(comm_, typeName, fName, data, fmt, ver, cmp, append);
    }
}


// ************************************************************************* //
