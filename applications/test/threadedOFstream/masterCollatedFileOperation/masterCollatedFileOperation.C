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

#include "masterCollatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "threadedCollatedOFstream.H"
#include "decomposedBlockData.H"
#include "registerSwitch.H"
#include "masterOFstream.H"
#include "OFstream.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(masterCollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        masterCollatedFileOperation,
        word
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::masterCollatedFileOperation::masterCollatedFileOperation
(
    const bool verbose
)
:
    collatedFileOperation(false)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName
            << " (maxThreadFileBufferSize " << maxThreadFileBufferSize
            << ')' << endl;

        if (maxThreadFileBufferSize == 0)
        {
            Info<< "         Threading not activated "
                   "since maxThreadFileBufferSize = 0." << nl
                << "         Writing may run slowly for large file sizes."
                << endl;
        }
        else
        {
            Info<< "         Threading activated "
                   "since maxThreadFileBufferSize > 0." << nl
                << "         Requires basic thread support"
                << " (MPI_THREAD_FUNNELED) enabled in MPI, "
                   "otherwise the simulation" << nl
                << "         may \"hang\".  If thread support cannot be "
                   "enabled, deactivate threading" << nl
                << "         by setting maxThreadFileBufferSize to 0 in "
                   "$FOAM_ETC/controlDict"
                << endl;
        }

        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::inotifyMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify" << endl;
        }

        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::timeStampMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::masterCollatedFileOperation::
~masterCollatedFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperations::masterCollatedFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    const Time& tm = io.time();
    const fileName& inst = io.instance();

    if (inst.isAbsolute() || !tm.processorCase())
    {
        mkDir(io.path());
        fileName pathName(io.objectPath());

        if (debug)
        {
            Pout<< "writeObject:"
                << " : For object : " << io.name()
                << " falling back to master-only output to " << io.path()
                << endl;
        }

        masterOFstream os
        (
            pathName,
            fmt,
            ver,
            cmp,
            false,
            valid
        );

        // If any of these fail, return (leave error handling to Ostream class)
        if (!os.good())
        {
            return false;
        }
        if (!io.writeHeader(os))
        {
            return false;
        }
        // Write the data to the Ostream
        if (!io.writeData(os))
        {
            return false;
        }
        IOobject::writeEndDivider(os);

        return true;
    }
    else
    {
        // Construct the equivalent processors/ directory
        fileName path(processorsPath(io, inst));

        mkDir(path);
        fileName pathName(path/io.name());

        if (io.global())
        {
            if (debug)
            {
                Pout<< "writeObject:" << " : For global object : " << io.name()
                    << " falling back to master-only output to " << pathName
                    << endl;
            }

            masterOFstream os
            (
                pathName,
                fmt,
                ver,
                cmp,
                false,
                valid
            );

            // If any of these fail, return (leave error handling to Ostream
            // class)
            if (!os.good())
            {
                return false;
            }
            if (!io.writeHeader(os))
            {
                return false;
            }
            // Write the data to the Ostream
            if (!io.writeData(os))
            {
                return false;
            }
            IOobject::writeEndDivider(os);

            return true;
        }
        else if (!Pstream::parRun())
        {
            // Special path for e.g. decomposePar. Append to
            // processors/ file
            if (debug)
            {
                Pout<< "writeObject:"
                    << " : For object : " << io.name()
                    << " appending to " << pathName << endl;
            }

            return appendObject(io, pathName, fmt, ver, cmp);
        }
        else
        {
            if (debug)
            {
                Pout<< "writeObject:"
                    << " : For object : " << io.name()
                    << " starting collating output to " << pathName << endl;
            }

            threadedCollatedOFstream os(writer_, pathName, fmt, ver, cmp);

            // If any of these fail, return (leave error handling to Ostream
            // class)
            if (!os.good())
            {
                return false;
            }
            if (Pstream::master() && !io.writeHeader(os))
            {
                return false;
            }
            // Write the data to the Ostream
            if (!io.writeData(os))
            {
                return false;
            }
            if (Pstream::master())
            {
                IOobject::writeEndDivider(os);
            }

            return true;
        }
    }
}


// ************************************************************************* //
