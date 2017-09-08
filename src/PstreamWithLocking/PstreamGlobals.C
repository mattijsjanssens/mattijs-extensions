/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "PstreamGlobals.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "OSspecific.H"
#include "error.H"
#include "IOstreams.H"
#include "OFstreamCollator.H"

#include <sys/syscall.h>
#include <mpi.h>
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Outstanding non-blocking operations.
//! \cond fileScope
DynamicList<MPI_Request> PstreamGlobals::outstandingRequests_;
//! \endcond

//// Max outstanding non-blocking operations.
////! \cond fileScope
//int PstreamGlobals::nRequests_ = 0;
////! \endcond

// Free'd non-blocking operations.
//! \cond fileScope
//DynamicList<label> PstreamGlobals::freedRequests_;
//! \endcond

// Max outstanding message tag operations.
//! \cond fileScope
int PstreamGlobals::nTags_ = 0;
//! \endcond

// Free'd message tags
//! \cond fileScope
DynamicList<int> PstreamGlobals::freedTags_;
//! \endcond


// Allocated communicators.
//! \cond fileScope
DynamicList<MPI_Comm> PstreamGlobals::MPICommunicators_;
DynamicList<MPI_Group> PstreamGlobals::MPIGroups_;
//! \endcond

void PstreamGlobals::checkCommunicator
(
    const label comm,
    const label otherProcNo
)
{
    if
    (
        comm < 0
     || comm >= PstreamGlobals::MPICommunicators_.size()
    )
    {
        FatalErrorInFunction
            << "otherProcNo:" << otherProcNo << " : illegal communicator "
            << comm << endl
            << "Communicator should be within range 0.."
            << PstreamGlobals::MPICommunicators_.size()-1 << abort(FatalError);
    }
}

void PstreamGlobals::checkThread(const label comm)
{
    if (comm == 1)
    {
        // Get current thread
        pid_t tid = syscall(SYS_gettid);

        if (OFstreamCollator::staticThread_ != tid)
        {
            Pout<< "** current:" << tid << endl;
            Pout<< "** stored:" << OFstreamCollator::staticThread_ << endl;
            Pout<< "** communicator:" << comm << endl;
            //error::printStack(Pout);
            //sleep(2000);
            //FatalErrorInFunction << "current:" << tid
            //    << " stored:" << commThread << " for communicator:" << comm
            //    << exit(FatalError);
        }
    }
}


// Thread lock
//! \cond fileScope
label PstreamGlobals::mutex_ = -1;
//! \endcond


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
