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
#include "masterOFstream.H"
#include "masterCollatingOFstream.H"
#include "decomposedBlockData.H"

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

Foam::fileName Foam::fileOperations::masterCollatingFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    pathType& searchType,
    word& newInstancePath
)
{
    if (io.name() == "controlDict" && io.instance() == "system")
    {
        // Special rules for system/controlDict since cannot use time() until
        // it has been loaded.

        fileName objectPath = io.path()/io.name();
        if (Foam::isFile(objectPath))
        {
            searchType = fileOperation::PROCESSORSOBJECT;
            return objectPath;
        }
        else if (Pstream::parRun())
        {
            // Search parent
            // Find out from case name whether a processor directory

            fileName globalCaseName;
            std::string::size_type pos = io.caseName().find("processor");
            if (pos != string::npos)
            {
                if (pos == 0)
                {
                    globalCaseName = ".";
                }
                else
                {
                    globalCaseName = io.caseName()(pos-1);
                }
            }

DebugVar(globalCaseName);

            if (globalCaseName.size())
            {
                fileName parentObjectPath =
                    io.rootPath()
                   /globalCaseName
                   /io.instance()
                   /io.db().dbDir()
                   /io.local()
                   /io.name();

DebugVar(parentObjectPath);

                if (Foam::isFile(parentObjectPath))
                {
                    searchType = fileOperation::PARENTOBJECT;
                    return parentObjectPath;
                }
            }
            return fileName::null;
        }
        else
        {
            return fileName::null;
        }
    }


DebugVar(io.objectPath());
DebugVar(checkGlobal);

    if (!io.time().processorCase())
    {
        return masterFileOperation::filePath
        (
            checkGlobal,
            io,
            searchType,
            newInstancePath
        );
    }


    newInstancePath = word::null;

    if (io.instance().isAbsolute())
    {
        fileName objectPath = io.instance()/io.name();
        if (Foam::isFile(objectPath))
        {
            searchType = fileOperation::ABSOLUTE;
            return objectPath;
        }
        else
        {
            searchType = fileOperation::NOTFOUND;
            return fileName::null;
        }
    }
    else
    {
        // 1. Check processors/
        fileName path = processorsPath(io, io.instance());
        fileName objectPath = path/io.name();
        if (Foam::isFile(objectPath))
        {
            searchType = fileOperation::PROCESSORSOBJECT;
            return objectPath;
        }
        else
        {
            // 2. Check local
            fileName localObjectPath = io.path()/io.name();

            if (Foam::isFile(localObjectPath))
            {
                searchType = fileOperation::OBJECT;
                return localObjectPath;
            }
        }



        // Any global checks
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
            fileName parentObjectPath =
                io.rootPath()/io.time().globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();

            if (Foam::isFile(parentObjectPath))
            {
                searchType = fileOperation::PARENTOBJECT;
                return parentObjectPath;
            }
        }


//        if (!Foam::isDir(path))
//        {
//            // 1. Check processors/
//            newInstancePath = io.time().findInstancePath
//            (
//                processorsCasePath(io),
//                instant(io.instance())
//            );
//
//            if (newInstancePath.size())
//            {
//                fileName fName(processorsPath(io, newInstancePath)/io.name());
//
//                if (Foam::isFile(fName))
//                {
//                    searchType = fileOperation::PROCESSORSFINDINSTANCE;
//                    return fName;
//                }
//            }
//        }
//
//        if (!Foam::isDir(io.path()))
//        {
//            // 2. Check local
//            newInstancePath = io.time().findInstancePath
//            (
//                instant(io.instance())
//            );
//
//            if (newInstancePath.size())
//            {
//                fileName fName
//                (
//                    io.rootPath()/io.caseName()
//                   /newInstancePath/io.db().dbDir()/io.local()/io.name()
//                );
//
//                if (Foam::isFile(fName))
//                {
//                    searchType = fileOperation::FINDINSTANCE;
//                    return fName;
//                }
//            }
//        }

        return fileName::null;
    }
}


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

Foam::fileName Foam::fileOperations::masterCollatingFileOperation::objectPath
(
    const IOobject& io
) const
{
    // Replacement for objectPath
    return masterFileOperation::objectPath
    (
        io,
        fileOperation::PROCESSORSOBJECT,
        io.instance()
    );
}


Foam::fileName Foam::fileOperations::masterCollatingFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io
) const
{
    fileName objPath;
    pathType searchType = fileOperation::NOTFOUND;
    word newInstancePath;
    if (Pstream::master())
    {
        objPath = filePath(checkGlobal, io, searchType, newInstancePath);
    }
    label masterType(searchType);
    Pstream::scatter(masterType);
    searchType = pathType(masterType);

    if
    (
        searchType == fileOperation::FINDINSTANCE
     || searchType == fileOperation::PROCESSORSFINDINSTANCE
    )
    {
        Pstream::scatter(newInstancePath);
    }

    if (!Pstream::master())
    {
        objPath = masterFileOperation::objectPath
        (
            io,
            searchType,
            newInstancePath
        );
    }
    return objPath;
}


Foam::autoPtr<Foam::Istream>
Foam::fileOperations::masterCollatingFileOperation::readStream
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

    autoPtr<Istream> isPtr;
    bool valid = true;
    if (UPstream::master())
    {
        isPtr.reset(new IFstream(fName));

        // Read header data (on copy)
        IOobject headerIO(io);
        headerIO.readHeader(isPtr());
        if (headerIO.headerClassName() != decomposedBlockData::typeName)
        {
            valid = false;
            isPtr.clear();
        }
    }

    Pstream::scatter(valid);

    if (!valid)
    {
        // Fall back

        if (IFstream::debug)
        {
            Pout<< "masterCollatingFileOperation::readStream :"
                << " for object : " << io.name()
                << " falling back to master-only reading from " << fName
                << endl;
        }

        return masterFileOperation::readStream(io, fName);
    }



    if (IFstream::debug)
    {
        Pout<< "masterCollatingFileOperation::writeObject :"
            << " for object : " << io.name()
            << " starting collating input from " << fName << endl;
    }

    // Read my data
    List<char> data;
    decomposedBlockData::readBlocks(isPtr, data, UPstream::scheduled);
    // TBD: remove extra copying
    string buf(data.begin(), data.size());
    autoPtr<Istream> realIsPtr(new IStringStream(buf));

    // Read header
    if (!io.readHeader(realIsPtr()))
    {
        FatalIOErrorInFunction(realIsPtr())
            << "problem while reading header for object " << io.name()
            << exit(FatalIOError);
    }

    return realIsPtr;
}


bool Foam::fileOperations::masterCollatingFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    autoPtr<Ostream> osPtr;
    if
    (
        io.global()
     || io.instance().isAbsolute()
     || !Pstream::parRun()
     || !io.time().processorCase()
    )
    {
        mkDir(io.path());
        fileName pathName(io.objectPath());

        if (IFstream::debug)
        {
            Pout<< "masterCollatingFileOperation::writeObject :"
                << " for object : " << io.name()
                << " falling back to master-only output to " << io.path()
                << endl;
        }
        osPtr.reset(new masterOFstream(pathName, fmt, ver, cmp));
    }
    else
    {
        // Construct the equivalent processors/ directory
        fileName path(processorsPath(io, io.instance()));

        if (IFstream::debug)
        {
            Pout<< "masterCollatingFileOperation::writeObject :"
                << " for object : " << io.name()
                << " starting collating output to " << path << endl;
        }

        mkDir(path);
        fileName pathName(path/io.name());

        osPtr.reset
        (
            new masterCollatingOFstream
            (
                pathName,
                UPstream::scheduled,
                fmt,
                ver
                //cmp
            )
        );
    }
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


// ************************************************************************* //
