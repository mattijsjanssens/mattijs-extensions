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

#include "masterFileServer.H"
#include "masterOFstream.H"
#include "Pstream.H"
#include "Time.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileServers
{
    defineTypeNameAndDebug(masterFileServer, 0);
    addToRunTimeSelectionTable(fileServer, masterFileServer, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::fileServers::masterFileServer::filePath
(
    const IOobject& io,
    pathType& searchType,
    word& newInstancePath
)
{
    newInstancePath = word::null;

    if (io.instance().isAbsolute())
    {
        fileName objectPath = io.instance()/io.name();
        if (Foam::isFile(objectPath))
        {
            searchType = fileServer::ABSOLUTE;
            return objectPath;
        }
        else
        {
            searchType = fileServer::NOTFOUND;
            return fileName::null;
        }
    }
    else
    {
        fileName path = io.path();
        fileName objectPath = path/io.name();

        if (Foam::isFile(objectPath))
        {
            searchType = fileServer::OBJECT;
            return objectPath;
        }
        else
        {
            if
            (
                io.time().processorCase()
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
                    searchType = fileServer::PARENTOBJECT;
                    return parentObjectPath;
                }
            }

            if (!Foam::isDir(path))
            {
                newInstancePath = io.time().findInstancePath
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
                        searchType = fileServer::FINDINSTANCE;
                        return fName;
                    }
                }
            }
        }

        return fileName::null;
    }
}


Foam::fileName Foam::fileServers::masterFileServer::objectPath
(
    const IOobject& io,
    const pathType& searchType,
    const word& instancePath
)
{
    // Replacement for IOobject::objectPath()

    switch (searchType)
    {
        case fileServer::ABSOLUTE:
        {
            return io.instance()/io.name();
        }
        break;

        case fileServer::OBJECT:
        {
            return io.path()/io.name();
        }
        break;

        case fileServer::PARENTOBJECT:
        {
            return
                io.rootPath()/io.time().globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case fileServer::FINDINSTANCE:
        {
            return
                io.rootPath()/io.caseName()
               /instancePath/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case fileServer::NOTFOUND:
        {
            return fileName::null;
        }
        break;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileServers::masterFileServer::masterFileServer()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileServers::masterFileServer::~masterFileServer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileServers::masterFileServer::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return masterFileOperation<mode_t, mkDirOp>(dir, mkDirOp(mode));
}


bool Foam::fileServers::masterFileServer::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return masterFileOperation<mode_t, chModOp>(fName, chModOp(mode));
}


mode_t Foam::fileServers::masterFileServer::mode(const fileName& fName) const
{
    return masterFileOperation<mode_t, modeOp>(fName, modeOp());
}


Foam::fileName::Type Foam::fileServers::masterFileServer::type
(
    const fileName& fName
) const
{
    return fileName::Type(masterFileOperation<label, typeOp>(fName, typeOp()));
}


bool Foam::fileServers::masterFileServer::exists
(
    const fileName& fName,
    const bool checkGzip
) const
{
    return masterFileOperation<bool, existsOp>(fName, existsOp(checkGzip));
}


bool Foam::fileServers::masterFileServer::isDir(const fileName& fName) const
{
    return masterFileOperation<bool, isDirOp>(fName, isDirOp());
}


bool Foam::fileServers::masterFileServer::isFile
(
    const fileName& fName,
    const bool checkGzip
) const
{
    return masterFileOperation<bool, isFileOp>(fName, isFileOp());
}


off_t Foam::fileServers::masterFileServer::fileSize(const fileName& fName) const
{
    return masterFileOperation<off_t, fileSizeOp>(fName, fileSizeOp());
}


time_t Foam::fileServers::masterFileServer::lastModified
(
    const fileName& fName
) const
{
    return masterFileOperation<time_t, lastModifiedOp>(fName, lastModifiedOp());
}


double Foam::fileServers::masterFileServer::highResLastModified
(
    const fileName& fName
) const
{
    return masterFileOperation<double, lastModifiedHROp>
    (
        fName,
        lastModifiedHROp()
    );
}


bool Foam::fileServers::masterFileServer::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return masterFileOperation<bool, mvBakOp>(fName, mvBakOp(ext));
}


bool Foam::fileServers::masterFileServer::rm(const fileName& fName) const
{
    return masterFileOperation<bool, rmOp>(fName, rmOp());
}


bool Foam::fileServers::masterFileServer::rmDir(const fileName& dir) const
{
    return masterFileOperation<bool, rmDirOp>(dir, rmDirOp());
}


Foam::fileNameList Foam::fileServers::masterFileServer::readDir
(
    const fileName& dir,
    const fileName::Type type,
    const bool filtergz
) const
{
    return masterFileOperation<fileNameList, readDirOp>
    (
        dir,
        readDirOp(type, filtergz)
    );
}


bool Foam::fileServers::masterFileServer::cp
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterFileOperation<bool, cpOp>(src, dst, cpOp());
}


bool Foam::fileServers::masterFileServer::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterFileOperation<bool, lnOp>(src, dst, lnOp());
}


bool Foam::fileServers::masterFileServer::mv
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterFileOperation<bool, mvOp>(src, dst, mvOp());
}


Foam::fileName Foam::fileServers::masterFileServer::filePath
(
    const IOobject& io
) const
{
    fileName objPath;
    pathType searchType = fileServer::NOTFOUND;
    word newInstancePath;
    if (Pstream::master())
    {
        objPath = filePath(io, searchType, newInstancePath);
    }
    label masterType(searchType);
    Pstream::scatter(masterType);
    searchType = pathType(masterType);
    if (searchType == fileServer::FINDINSTANCE)
    {
        Pstream::scatter(newInstancePath);
    }

    if (!Pstream::master())
    {
        objPath = objectPath(io, searchType, newInstancePath);
    }
    return objPath;
}


Foam::autoPtr<Foam::Istream>
Foam::fileServers::masterFileServer::NewIFstream(const fileName& filePath) const
{
    if (Pstream::parRun())
    {
        // Insert logic of filePath. We assume that if a file is absolute
        // on the master it is absolute also on the slaves etc.

        List<fileName> filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = filePath;
        Pstream::gatherList(filePaths);

        PstreamBuffers pBufs(Pstream::nonBlocking);

        if (Pstream::master())
        {
            bool uniform = true;
            const fileName& object0 = filePaths[0];

            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
                if (filePaths[proci] != object0)
                {
                    uniform = false;
                    break;
                }
            }

            if (uniform)
            {
                if (IFstream::debug)
                {
                    Pout<< "Opening global file " << object0 << endl;
                }

                std::ifstream is(object0);
                // get length of file:
                is.seekg(0, ios_base::end);
                std::streamoff count = is.tellg();
                is.seekg(0, ios_base::beg);

                if (IFstream::debug)
                {
                    Pout<< "From " << object0
                        <<  " reading " << label(count) << " bytes" << endl;
                }
                List<char> buf(count);
                is.read(buf.begin(), count);

                for (label proci = 1; proci < Pstream::nProcs(); proci++)
                {
                    UOPstream os(proci, pBufs);
                    os.write(buf.begin(), count);
                }
            }
            else
            {
                for (label proci = 1; proci < Pstream::nProcs(); proci++)
                {
                    if (IFstream::debug)
                    {
                        Pout<< "For processor " << proci
                            << " opening " << filePaths[proci] << endl;
                    }

                    std::ifstream is(filePaths[proci]);
                    // Get length of file
                    is.seekg(0, ios_base::end);
                    std::streamoff count = is.tellg();
                    is.seekg(0, ios_base::beg);

                    if (IFstream::debug)
                    {
                        Pout<< "From " << filePaths[proci]
                            <<  " reading " << label(count) << " bytes" << endl;
                    }
                    List<char> buf(count);
                    is.read(buf.begin(), count);

                    UOPstream os(proci, pBufs);
                    os.write(buf.begin(), count);
                }
            }
        }


        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (Pstream::master())
        {
            // Read myself
            return autoPtr<Istream>
            (
                new IFstream(filePaths[Pstream::masterNo()])
            );
        }
        else
        {
            if (IFstream::debug)
            {
                Pout<< "Reading " << filePath
                    << " from processor " << Pstream::masterNo()
                    << endl;
            }

            UIPstream is(Pstream::masterNo(), pBufs);
            List<char> buf(recvSizes[Pstream::masterNo()]);
            is.read(buf.begin(), buf.size());

            if (IFstream::debug)
            {
                Pout<< "Done reading " << buf.size() << " bytes" << endl;
            }

            // Note: IPstream is not an IStream so use a IStringStream to
            //       convert the buffer.
            return autoPtr<Istream>(new IStringStream(buf.begin()));
        }
    }
    else
    {
        // Read myself
        return autoPtr<Istream>(new IFstream(filePath));
    }
}


Foam::autoPtr<Foam::Ostream> Foam::fileServers::masterFileServer::NewOFstream
(
    const fileName& pathName,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return autoPtr<Ostream>(new masterOFstream(pathName, fmt, ver, cmp));
}


// ************************************************************************* //
