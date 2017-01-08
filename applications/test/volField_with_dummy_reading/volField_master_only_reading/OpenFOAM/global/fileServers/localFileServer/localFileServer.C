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

#include "localFileServer.H"
#include "Time.H"
#include "IFstream.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileServers
{
    defineTypeNameAndDebug(localFileServer, 0);
    addToRunTimeSelectionTable(fileServer, localFileServer, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileServers::localFileServer::localFileServer()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileServers::localFileServer::~localFileServer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileServers::localFileServer::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return Foam::mkDir(dir, mode);
}


bool Foam::fileServers::localFileServer::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return Foam::chMod(fName, mode);
}


mode_t Foam::fileServers::localFileServer::mode(const fileName& fName) const
{
    return Foam::mode(fName);
}


Foam::fileName::Type Foam::fileServers::localFileServer::type
(
    const fileName& fName
) const
{
    return Foam::type(fName);
}


bool Foam::fileServers::localFileServer::exists
(
    const fileName& fName,
    const bool checkGzip
) const
{
    return Foam::exists(fName, checkGzip);
}


bool Foam::fileServers::localFileServer::isDir(const fileName& fName) const
{
    return Foam::isDir(fName);
}


bool Foam::fileServers::localFileServer::isFile
(
    const fileName& fName,
    const bool checkGzip
) const
{
    return Foam::isFile(fName, checkGzip);
}


off_t Foam::fileServers::localFileServer::fileSize(const fileName& fName) const
{
    return Foam::fileSize(fName);
}


time_t Foam::fileServers::localFileServer::lastModified
(
    const fileName& fName
) const
{
    return Foam::lastModified(fName);
}


double Foam::fileServers::localFileServer::highResLastModified
(
    const fileName& fName
) const
{
    return Foam::highResLastModified(fName);
}


bool Foam::fileServers::localFileServer::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return Foam::mvBak(fName, ext);
}


bool Foam::fileServers::localFileServer::rm(const fileName& fName) const
{
    return Foam::rm(fName);
}


bool Foam::fileServers::localFileServer::rmDir(const fileName& dir) const
{
    return Foam::rmDir(dir);
}


Foam::fileNameList Foam::fileServers::localFileServer::readDir
(
    const fileName& dir,
    const fileName::Type type,
    const bool filtergz
) const
{
    return Foam::readDir(dir, type, filtergz);
}


bool Foam::fileServers::localFileServer::cp
(
    const fileName& src,
    const fileName& dst
) const
{
    return Foam::cp(src, dst);
}


bool Foam::fileServers::localFileServer::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return Foam::ln(src, dst);
}


bool Foam::fileServers::localFileServer::mv
(
    const fileName& src,
    const fileName& dst
) const
{
    return Foam::mv(src, dst);
}


Foam::fileName Foam::fileServers::localFileServer::filePath
(
    const IOobject& io
) const
{
    return io.filePath();
}


Foam::autoPtr<Foam::Istream> Foam::fileServers::localFileServer::objectStream
(
    const fileName& fName
) const
{
    if (fName.size())
    {
        autoPtr<Istream> isPtr = NewIFstream(fName);

        if (isPtr->good())
        {
            return isPtr;
        }
    }
    return autoPtr<Istream>(nullptr);
}


Foam::autoPtr<Foam::Istream>
Foam::fileServers::localFileServer::NewIFstream(const fileName& filePath) const
{
    return autoPtr<Istream>(new IFstream(filePath));
}


Foam::autoPtr<Foam::Ostream> Foam::fileServers::localFileServer::NewOFstream
(
    const fileName& pathName,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return autoPtr<Ostream>(new OFstream(pathName, fmt, ver, cmp));
}


// ************************************************************************* //
