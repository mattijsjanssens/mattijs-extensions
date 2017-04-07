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

#include "collatedFileOperation2.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "OFstreamWriter.H"
#include "threadedOFstream.H"
#include "masterOFstream.H"
#include "masterCollatingOFstream.H"
#include "decomposedBlockData.H"
#include "OFstream.H"
#include "registerSwitch.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(collatedFileOperation2, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        collatedFileOperation2,
        word
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation2::collatedFileOperation2
(
    const bool verbose
)
:
    collatedFileOperation(false)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName
            << " (maxBufferSize " << maxBufferSize
            << "  maxThreadBufferSize " << maxThreadBufferSize
            << ')' << endl;
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

Foam::fileOperations::collatedFileOperation2::~collatedFileOperation2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperations::collatedFileOperation2::writeObject
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
        return collatedFileOperation2::writeObject(io, fmt, ver, cmp, valid);
    }
    else
    {
        // Construct the equivalent processors/ directory
        fileName path(processorsPath(io, inst));

        mkDir(path);
        fileName pathName(path/io.name());

        if (io.global() || !Pstream::parRun()))
        {
            return collatedFileOperation2::writeObject(io, fmt, ver, cmp, valid);
        }
        else
        {
            if (debug)
            {
                Pout<< "writeObject:"
                    << " : For object : " << io.name()
                    << " starting collating output to " << path << endl;
            }
//XXXXXXX
            osPtr.reset
            (
                new masterCollatingOFstream
                (
                    (
                        off_t(maxThreadBufferSize) > 0
                     ? &writeServer_
                     : nullptr
                    ),
                    pathName,
                    fmt,
                    ver,
                    cmp
                )
            );

            // Suppress header on slaves
            writeHeader = Pstream::master();

            Ostream& os = osPtr();

            // If any of these fail, return (leave error handling to Ostream class)
            if (!os.good())
            {
                return false;
            }

            if (writeHeader && !io.writeHeader(os))
            {
                return false;
            }

            // Write the data to the Ostream
            if (!io.writeData(os))
            {
                return false;
            }

            if (writeHeader)
            {
                IOobject::writeEndDivider(os);
            }

    return true;
}


// ************************************************************************* //
