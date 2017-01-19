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

#include "localFileOperation.H"
#include "Time.H"
#include "IFstream.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(localFileOperation, 0);
    addToRunTimeSelectionTable(fileOperation, localFileOperation, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::localFileOperation::localFileOperation()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::localFileOperation::~localFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperations::localFileOperation::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return Foam::mkDir(dir, mode);
}


bool Foam::fileOperations::localFileOperation::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return Foam::chMod(fName, mode);
}


mode_t Foam::fileOperations::localFileOperation::mode
(
    const fileName& fName
) const
{
    return Foam::mode(fName);
}


Foam::fileName::Type Foam::fileOperations::localFileOperation::type
(
    const fileName& fName
) const
{
    return Foam::type(fName);
}


bool Foam::fileOperations::localFileOperation::exists
(
    const fileName& fName,
    const bool checkGzip
) const
{
    return Foam::exists(fName, checkGzip);
}


bool Foam::fileOperations::localFileOperation::isDir
(
    const fileName& fName
) const
{
    return Foam::isDir(fName);
}


bool Foam::fileOperations::localFileOperation::isFile
(
    const fileName& fName,
    const bool checkGzip
) const
{
    return Foam::isFile(fName, checkGzip);
}


off_t Foam::fileOperations::localFileOperation::fileSize
(
    const fileName& fName
) const
{
    return Foam::fileSize(fName);
}


time_t Foam::fileOperations::localFileOperation::lastModified
(
    const fileName& fName
) const
{
    return Foam::lastModified(fName);
}


double Foam::fileOperations::localFileOperation::highResLastModified
(
    const fileName& fName
) const
{
    return Foam::highResLastModified(fName);
}


bool Foam::fileOperations::localFileOperation::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return Foam::mvBak(fName, ext);
}


bool Foam::fileOperations::localFileOperation::rm(const fileName& fName) const
{
    return Foam::rm(fName);
}


bool Foam::fileOperations::localFileOperation::rmDir(const fileName& dir) const
{
    return Foam::rmDir(dir);
}


Foam::fileNameList Foam::fileOperations::localFileOperation::readDir
(
    const fileName& dir,
    const fileName::Type type,
    const bool filtergz
) const
{
    return Foam::readDir(dir, type, filtergz);
}


bool Foam::fileOperations::localFileOperation::cp
(
    const fileName& src,
    const fileName& dst
) const
{
    return Foam::cp(src, dst);
}


bool Foam::fileOperations::localFileOperation::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return Foam::ln(src, dst);
}


bool Foam::fileOperations::localFileOperation::mv
(
    const fileName& src,
    const fileName& dst
) const
{
    return Foam::mv(src, dst);
}


Foam::fileName Foam::fileOperations::localFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io
) const
{
    if (io.instance().isAbsolute())
    {
        fileName objectPath = io.instance()/io.name();

        if (Foam::isFile(objectPath))
        {
            return objectPath;
        }
        else
        {
            return fileName::null;
        }
    }
    else
    {
        fileName path = io.path();
        fileName objectPath = path/io.name();

        if (Foam::isFile(objectPath))
        {
            return objectPath;
        }
        else
        {
            if
            (
                checkGlobal
             && io.time().processorCase()
             && (
                    io.instance() == io.time().system()
                 || io.instance() == io.time().constant()
                )
            )
            {
                // Constant & system can come from global case

                fileName parentObjectPath =
                    io.rootPath()/io.time().globalCaseName()
                   /io.instance()/io.db().dbDir()/io.local()/io.name();

                if (Foam::isFile(parentObjectPath))
                {
                    return parentObjectPath;
                }
            }

            // Check for approximately same time
            if (!Foam::isDir(path))
            {
                word newInstancePath = io.time().findInstancePath
                (
                    instant(io.instance())
                );

                if (newInstancePath.size())
                {
                    fileName fName
                    (
                        io.rootPath()/io.caseName()
                       /newInstancePath/io.db().dbDir()/io.local()/io.name()
                    );

                    if (Foam::isFile(fName))
                    {
                        return fName;
                    }
                }
            }
        }

        return fileName::null;
    }
}


Foam::autoPtr<Foam::Istream>
Foam::fileOperations::localFileOperation::objectStream
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
Foam::fileOperations::localFileOperation::readStream
(
    regIOobject& io,
    const fileName& fName
) const
{
    if (!fName.size())
    {
        FatalErrorInFunction
            << "empty file name" << exit(FatalError);
    }

    autoPtr<Istream> isPtr = NewIFstream(fName);

    if (!isPtr.valid() || !isPtr->good())
    {
        FatalIOError
        (
            "localFileOperation::readStream()",
            __FILE__,
            __LINE__,
            fName,
            0
        )   << "cannot open file"
            << exit(FatalIOError);
    }
    else if (!io.readHeader(isPtr()))
    {
        FatalIOErrorInFunction(isPtr())
            << "problem while reading header for object " << io.name()
            << exit(FatalIOError);
    }

    return isPtr;
}


Foam::autoPtr<Foam::Istream>
Foam::fileOperations::localFileOperation::NewIFstream
(
    const fileName& filePath
) const
{
    return autoPtr<Istream>(new IFstream(filePath));
}


Foam::autoPtr<Foam::Ostream>
Foam::fileOperations::localFileOperation::NewOFstream
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
