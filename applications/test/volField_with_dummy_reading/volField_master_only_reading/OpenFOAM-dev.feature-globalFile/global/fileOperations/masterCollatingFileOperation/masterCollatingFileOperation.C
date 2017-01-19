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

#include "masterCollatingFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "IFstream.H"
#include "masterCollatingOFstream.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(masterCollatingFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        masterCollatingFileOperation,
        word
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::masterCollatingFileOperation::
masterCollatingFileOperation()
:
    masterFileOperation()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::masterCollatingFileOperation::
~masterCollatingFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//Foam::autoPtr<Foam::Istream>
//Foam::fileOperations::masterFileOperation::readStream
//(
//    regIOobject& io,
//    const fileName& fName
//) const
//{
//    if (!fName.size())
//    {
//        FatalErrorInFunction
//            << "empty file name" << exit(FatalError);
//    }
//
//    autoPtr<Istream> isPtr;
//    if (UPstream::master())
//    {
//        isPtr.reset(new IFstream(fName));
//        io.readHeader(isPtr());
//    }
//    return readBlocks(isPtr, *this, commsType_);
//
//    IOobject
//}


bool Foam::fileOperations::masterCollatingFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    mkDir(io.path());
    fileName pathName(io.objectPath());

    autoPtr<Ostream> osPtr(NewOFstream(pathName, fmt, ver, cmp));
    Ostream& os = osPtr();

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


Foam::autoPtr<Foam::Ostream>
Foam::fileOperations::masterCollatingFileOperation::NewOFstream
(
    const fileName& pathName,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return autoPtr<Ostream>
    (
        new masterCollatingOFstream
        (
            pathName,
            UPstream::scheduled,
            fmt,
            ver,
            cmp
        )
    );
}


// ************************************************************************* //
