/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "allReduce.H"
#include <fstream>
#include <cstdlib>
#include <cctype>

#include <stdio.h>
#include <unistd.h>
#include <dirent.h>
#include <pwd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <netdb.h>
#include <dlfcn.h>
#include <link.h>
#include <sys/syscall.h>

#include <netinet/in.h>
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, class BinaryOp>
void Foam::allReduce
(
    Type& Value,
    int MPICount,
    MPI_Datatype MPIType,
    MPI_Op MPIOp,
    const BinaryOp& bop,
    const int tag,
    const label communicator
)
{

if (Pstream::debug)
{
//pid_t tid = gettid();
pid_t tid = syscall(SYS_gettid);
//pthread_id_np_t   tid;
//tid = pthread_getthreadid_np();
std::cout<< "My thread:" << tid << " mutex:" << PstreamGlobals::mutex_
    << " on comm:" << communicator << std::endl;
}

    if (!UPstream::parRun())
    {
        return;
    }

    PstreamGlobals::checkThread(communicator);
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "UIPstream::allReduce :"
            << " tag:" << tag << " comm:" << communicator
            << Foam::endl;
        error::printStack(Pout);
    }

    PstreamGlobals::timer_.cpuTimeIncrement();
//     struct timespec t;
//     clock_gettime(CLOCK_REALTIME, &t);

    if (UPstream::nProcs(communicator) <= UPstream::nProcsSimpleSum)
    {
        if (UPstream::master(communicator))
        {
            for
            (
                int slave=UPstream::firstSlave();
                slave<=UPstream::lastSlave(communicator);
                slave++
            )
            {
                Type value;

                //lockMutex(PstreamGlobals::mutex_);
                if
                (
                    MPI_Recv
                    (
                        &value,
                        MPICount,
                        MPIType,
                        slave,  //UPstream::procID(slave),
                        tag,
                        PstreamGlobals::MPICommunicators_[communicator],
                        MPI_STATUS_IGNORE
                    )
                )
                {
                    FatalErrorInFunction
                        << "MPI_Recv failed"
                        << Foam::abort(FatalError);
                }
                //unlockMutex(PstreamGlobals::mutex_);

                Value = bop(Value, value);
            }
        }
        else
        {
            //lockMutex(PstreamGlobals::mutex_);
            if
            (
                MPI_Send
                (
                    &Value,
                    MPICount,
                    MPIType,
                    UPstream::masterNo(),//UPstream::procID(masterNo()),
                    tag,
                    PstreamGlobals::MPICommunicators_[communicator]
                )
            )
            {
                FatalErrorInFunction
                    << "MPI_Send failed"
                    << Foam::abort(FatalError);
            }
            //unlockMutex(PstreamGlobals::mutex_);
        }


        if (UPstream::master(communicator))
        {
            for
            (
                int slave=UPstream::firstSlave();
                slave<=UPstream::lastSlave(communicator);
                slave++
            )
            {
                //lockMutex(PstreamGlobals::mutex_);
                if
                (
                    MPI_Send
                    (
                        &Value,
                        MPICount,
                        MPIType,
                        slave,      //UPstream::procID(slave),
                        tag,
                        PstreamGlobals::MPICommunicators_[communicator]
                    )
                )
                {
                    FatalErrorInFunction
                        << "MPI_Send failed"
                        << Foam::abort(FatalError);
                }
                //unlockMutex(PstreamGlobals::mutex_);
            }
        }
        else
        {
            //lockMutex(PstreamGlobals::mutex_);
            if
            (
                MPI_Recv
                (
                    &Value,
                    MPICount,
                    MPIType,
                    UPstream::masterNo(),//UPstream::procID(masterNo()),
                    tag,
                    PstreamGlobals::MPICommunicators_[communicator],
                    MPI_STATUS_IGNORE
                )
            )
            {
                FatalErrorInFunction
                    << "MPI_Recv failed"
                    << Foam::abort(FatalError);
            }
            //unlockMutex(PstreamGlobals::mutex_);
        }
    }
    else
    {
        //lockMutex(PstreamGlobals::mutex_);
        Type sum;
        MPI_Allreduce
        (
            &Value,
            &sum,
            MPICount,
            MPIType,
            MPIOp,
            PstreamGlobals::MPICommunicators_[communicator]
        );
        //unlockMutex(PstreamGlobals::mutex_);
        Value = sum;
    }

    PstreamGlobals::reduceTime_ += PstreamGlobals::timer_.cpuTimeIncrement();
}


// ************************************************************************* //
