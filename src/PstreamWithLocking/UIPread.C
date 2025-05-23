/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Description
    Read from UIPstream

\*---------------------------------------------------------------------------*/

#include "UIPstream.H"
#include "PstreamGlobals.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "SubList.H"
#include "allReduce.H"

#include <sys/syscall.h>
#include "OSspecific.H"
#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <csignal>

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UIPstream::UIPstream
(
    const commsTypes commsType,
    const int fromProcNo,
    DynamicList<char>& externalBuf,
    label& externalBufPosition,
    const int tag,
    const label comm,
    const bool clearAtEnd,
    streamFormat format,
    versionNumber version
)
:
    UPstream(commsType),
    Istream(format, version),
    fromProcNo_(fromProcNo),
    externalBuf_(externalBuf),
    externalBufPosition_(externalBufPosition),
    tag_(tag),
    comm_(comm),
    clearAtEnd_(clearAtEnd),
    messageSize_(0)
{
    PstreamGlobals::checkThread(comm);

    setOpened();
    setGood();

    if (commsType == commsTypes::nonBlocking)
    {
        // Message is already received into externalBuf
    }
    else
    {
        MPI_Status status;

        label wantedSize = externalBuf_.capacity();

        if (debug)
        {
            pid_t tid = syscall(SYS_gettid);
            Pout<< "UIPstream::UIPstream : read from:" << fromProcNo
                << " tag:" << tag << " comm:" << comm_
                << " thread:" << tid
                << " wanted size:" << wantedSize
                << Foam::endl;
        }


        // If the buffer size is not specified, probe the incomming message
        // and set it
        if (!wantedSize)
        {
            //lockMutex(PstreamGlobals::mutex_);
            MPI_Probe
            (
                fromProcNo_,
                tag_,
                PstreamGlobals::MPICommunicators_[comm_],
                &status
            );
            MPI_Get_count(&status, MPI_BYTE, &messageSize_);
            //unlockMutex(PstreamGlobals::mutex_);

            externalBuf_.setCapacity(messageSize_);
            wantedSize = messageSize_;

            if (debug)
            {
                Pout<< "UIPstream::UIPstream : probed size:" << wantedSize
                    << Foam::endl;
            }
        }

        messageSize_ = UIPstream::read
        (
            commsType,
            fromProcNo_,
            externalBuf_.begin(),
            wantedSize,
            tag_,
            comm_
        );

        // Set addressed size. Leave actual allocated memory intact.
        externalBuf_.setSize(messageSize_);

        if (!messageSize_)
        {
            setEof();
        }
    }
}


Foam::UIPstream::UIPstream(const int fromProcNo, PstreamBuffers& buffers)
:
    UPstream(buffers.commsType_),
    Istream(buffers.format_, buffers.version_),
    fromProcNo_(fromProcNo),
    externalBuf_(buffers.recvBuf_[fromProcNo]),
    externalBufPosition_(buffers.recvBufPos_[fromProcNo]),
    tag_(buffers.tag_),
    comm_(buffers.comm_),
    clearAtEnd_(true),
    messageSize_(0)
{
    if
    (
        commsType() != UPstream::commsTypes::scheduled
     && !buffers.finishedSendsCalled_
    )
    {
        FatalErrorInFunction
            << "PstreamBuffers::finishedSends() never called." << endl
            << "Please call PstreamBuffers::finishedSends() after doing"
            << " all your sends (using UOPstream) and before doing any"
            << " receives (using UIPstream)" << Foam::exit(FatalError);
    }

    PstreamGlobals::checkThread(comm_);

    setOpened();
    setGood();

    if (commsType() == commsTypes::nonBlocking)
    {
        // Message is already received into externalBuf
        messageSize_ = buffers.recvBuf_[fromProcNo].size();

        if (debug)
        {
            pid_t tid = syscall(SYS_gettid);
            Pout<< "UIPstream::UIPstream PstreamBuffers :"
                << " fromProcNo:" << fromProcNo
                << " tag:" << tag_ << " comm:" << comm_
                << " thread:" << tid
                << " receive buffer size:" << messageSize_
                << Foam::endl;
        }
    }
    else
    {
        MPI_Status status;

        label wantedSize = externalBuf_.capacity();

        if (debug)
        {
            pid_t tid = syscall(SYS_gettid);
            Pout<< "UIPstream::UIPstream PstreamBuffers :"
                << " read from:" << fromProcNo
                << " tag:" << tag_ << " comm:" << comm_
                << " thread:" << tid
                << " wanted size:" << wantedSize
                << Foam::endl;
        }

        // If the buffer size is not specified, probe the incomming message
        // and set it
        if (!wantedSize)
        {
            PstreamGlobals::timer_.cpuTimeIncrement();
            //lockMutex(PstreamGlobals::mutex_);
            MPI_Probe
            (
                fromProcNo_,
                tag_,
                PstreamGlobals::MPICommunicators_[comm_],
                &status
            );
            MPI_Get_count(&status, MPI_BYTE, &messageSize_);
            //unlockMutex(PstreamGlobals::mutex_);
            PstreamGlobals::waitTime_ +=
                PstreamGlobals::timer_.cpuTimeIncrement();

            externalBuf_.setCapacity(messageSize_);
            wantedSize = messageSize_;

            if (debug)
            {
                Pout<< "UIPstream::UIPstream PstreamBuffers : probed size:"
                    << wantedSize << Foam::endl;
            }
        }

        messageSize_ = UIPstream::read
        (
            commsType(),
            fromProcNo_,
            externalBuf_.begin(),
            wantedSize,
            tag_,
            comm_
        );

        // Set addressed size. Leave actual allocated memory intact.
        externalBuf_.setSize(messageSize_);

        if (!messageSize_)
        {
            setEof();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UIPstream::read
(
    const commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize,
    const int tag,
    const label communicator
)
{
    if (debug)
    {
        pid_t tid = syscall(SYS_gettid);
        Pout<< "UIPstream::read : starting read from:" << fromProcNo
            << " tag:" << tag << " comm:" << communicator
            << " thread:" << tid
            << " wanted size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << Foam::endl;
    }
    PstreamGlobals::checkThread(communicator);
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "UIPstream::read : starting read from:" << fromProcNo
            << " tag:" << tag << " comm:" << communicator
            << " wanted size:" << label(bufSize)
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Pout);
    }

    if (commsType == commsTypes::blocking || commsType == commsTypes::scheduled)
    {
        MPI_Status status;

        PstreamGlobals::timer_.cpuTimeIncrement();
        //lockMutex(PstreamGlobals::mutex_);
        if
        (
            MPI_Recv
            (
                buf,
                bufSize,
                MPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &status
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Recv cannot receive incomming message"
                << Foam::abort(FatalError);

            return 0;
        }
        PstreamGlobals::waitTime_ +=
            PstreamGlobals::timer_.cpuTimeIncrement();


        // Check size of message read

        int messageSize;
        MPI_Get_count(&status, MPI_BYTE, &messageSize);

        //unlockMutex(PstreamGlobals::mutex_);

        if (debug)
        {
            pid_t tid = syscall(SYS_gettid);
            Pout<< "UIPstream::read : finished read from:" << fromProcNo
                << " tag:" << tag << " read size:" << label(bufSize)
                << " thread:" << tid
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << Foam::endl;
        }

        if (messageSize > bufSize)
        {
            FatalErrorInFunction
                << "buffer (" << label(bufSize)
                << ") not large enough for incomming message ("
                << messageSize << ')'
                << Foam::abort(FatalError);
        }

        return messageSize;
    }
    else if (commsType == commsTypes::nonBlocking)
    {
        MPI_Request request;

        //lockMutex(PstreamGlobals::mutex_);
        if
        (
            MPI_Irecv
            (
                buf,
                bufSize,
                MPI_BYTE,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Recv cannot start non-blocking receive"
                << Foam::abort(FatalError);

            return 0;
        }
        //unlockMutex(PstreamGlobals::mutex_);

        if (debug)
        {
            pid_t tid = syscall(SYS_gettid);
            Pout<< "UIPstream::read : started read from:" << fromProcNo
                << " tag:" << tag << " read size:" << label(bufSize)
                << " comm:" << communicator
                << " commsType:" << UPstream::commsTypeNames[commsType]
                << " thread:" << tid
                << " request:" << PstreamGlobals::outstandingRequests_.size()
                << Foam::endl;
        }

        PstreamGlobals::outstandingRequests_.append(request);

        // Assume the message is completely received.
        return bufSize;
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << int(commsType)
            << Foam::abort(FatalError);

        return 0;
    }
}


// ************************************************************************* //
