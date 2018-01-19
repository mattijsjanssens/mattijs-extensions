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
// #include "Time.H"
// #include "threadedCollatedOFstream.H"
// #include "decomposedBlockData.H"
// #include "registerSwitch.H"
// #include "masterOFstream.H"
// #include "OFstream.H"
// #include "collatedFileOperation.H"

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

// Foam::labelPair Foam::fileOperations::mixedCollatedFileOperation::commsGroup
// (
//     const label nProcs,
//     const label proci
// )
// {
//     label n = 2;    //1;
//     label master = (proci / n) * n;
//     return labelPair(master, n);
// }


Foam::labelList Foam::fileOperations::mixedCollatedFileOperation::subRanks
(
    const label n
)
{
    const string myHostName(hostName());

    stringList hosts(Pstream::nProcs());
    hosts[Pstream::myProcNo()] = myHostName;
    Pstream::gatherList(hosts);
    Pstream::scatterList(hosts);

    // Collect procs with same hostname
    DynamicList<label> subRanks(64);
    forAll(hosts, proci)
    {
        if (hosts[proci] == myHostName)
        {
            subRanks.append(proci);
        }
    }
    return subRanks;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::mixedCollatedFileOperation::mixedCollatedFileOperation
(
    const bool verbose
)
:
    collatedFileOperation
    (
        UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            subRanks(Pstream::nProcs())
        ),
        true
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::mixedCollatedFileOperation::~mixedCollatedFileOperation()
{
    if (comm_ != -1)
    {
        UPstream::freeCommunicator(comm_);
    }
}


// ************************************************************************* //
