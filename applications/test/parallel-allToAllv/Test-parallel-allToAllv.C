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

Application
    Test-parallel-nonBlocking

Description
    Test for various non-blocking parallel routines.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "mapDistribute.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
// #include "vector.H"
// #include "IOstreams.H"
// #include "Random.H"
// #include "Tuple2.H"
// #include "PstreamBuffers.H"

#define MPICH_SKIP_MPICXX
#include "/home/penfold2/mattijs/OpenFOAM/ThirdParty-dev/platforms/linux64Gcc/mpich-3.2/include/mpi.h"
#include "UPstream.H"
#include "/home/penfold2/mattijs/OpenFOAM/OpenFOAM-dev/src/Pstream/mpi/PstreamGlobals.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void allToAll
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,

    const label communicator
)
{
    label np = UPstream::nProcs(communicator);

    if
    (
        sendSizes.size() != np
     || sendOffsets.size() != np
     || recvSizes.size() != np
     || recvOffsets.size() != np
    )
    {
        FatalErrorInFunction
            << "Size of sendSize " << sendSizes.size()
            << ", sendOffsets " << sendOffsets.size()
            << ", recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if
    (
        MPI_Alltoallv
        (
            sendData,
            sendSizes.begin(),
            sendOffsets.begin(),
            MPI_BYTE,
            recvData,
            recvSizes.begin(),
            recvOffsets.begin(),
            MPI_BYTE,
            PstreamGlobals::MPICommunicators_[communicator]
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Alltoallv failed for sendSizes " << sendSizes
            << " recvSizes " << recvSizes
            << " communicator " << communicator
            << Foam::abort(FatalError);
    }
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "End\n" << endl;

    List<int> sendSizes(Pstream::nProcs());
    List<int> sendOffsets(Pstream::nProcs());

    List<int> recvSizes(Pstream::nProcs());
    List<int> recvOffsets(Pstream::nProcs());

    // E.g. send my pids over to all other processors

    labelList allPids(Pstream::nProcs());
    allPids[Pstream::myProcNo()] = pid();

    DebugVar(allPids[Pstream::myProcNo()]);

    // Replacement for gatherList
    {
        char* data = reinterpret_cast<char*>(allPids.begin());

        forAll(sendSizes, proci)
        {
            if (proci == Pstream::myProcNo())
            {
                sendSizes[proci] = 0;
                sendOffsets[proci] = 0;

                recvSizes[proci] = 0;
                recvOffsets[proci] = 0;
            }
            else
            {
                sendSizes[proci] = sizeof(label);
                sendOffsets[proci] =
                    reinterpret_cast<char*>(&allPids[Pstream::myProcNo()])
                  - data;

                recvSizes[proci] = sizeof(label);
                recvOffsets[proci] =
                    reinterpret_cast<char*>(&allPids[proci])
                  - data;
            }
        }

        allToAll
        (
            data,
            sendSizes,
            sendOffsets,

            data,
            recvSizes,
            recvOffsets,

            Pstream::worldComm
        );
    }
    DebugVar(allPids);



    return 0;
}


// ************************************************************************* //
