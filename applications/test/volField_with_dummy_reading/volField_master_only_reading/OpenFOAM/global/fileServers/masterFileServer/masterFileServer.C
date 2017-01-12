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
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "IFstream.H"
#include "masterOFstream.H"
#include "masterCollatingOFstream.H"

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


bool Foam::fileServers::masterFileServer::uniformFile
(
    const fileNameList& filePaths
)
{
    const fileName& object0 = filePaths[0];

    for (label i = 1; i < filePaths.size(); i++)
    {
        if (filePaths[i] != object0)
        {
            return false;
        }
    }
    return true;
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


//Foam::autoPtr<Foam::Istream> Foam::fileServers::masterFileServer::objectStream
//(
//    const fileName& fName
//) const
//{
//    if (fName.size())
//    {
//        autoPtr<Istream> isPtr = NewIFstream(fName);
//
//        if (isPtr->good())
//        {
//            return isPtr;
//        }
//    }
//    return autoPtr<Istream>(nullptr);
//}


//Foam::autoPtr<Foam::Istream> Foam::fileServers::masterFileServer::readStream
//(
//    regIOobject& io,
//    const fileName& fName
//) const
//{
//    //autoPtr<Istream> isPtr = objectStream(fName);
//
//    if (!fName.size())
//    {
//        FatalErrorInFunction
//            << "empty file name" << exit(FatalError);
//    }
//
//    autoPtr<Istream> isPtr = NewIFstream(fName);
//
//    if (!isPtr.valid() || !isPtr->good())
//    {
//        FatalIOError
//        (
//            "masterFileServer::readStream()",
//            __FILE__,
//            __LINE__,
//            fName,
//            0
//        )   << "cannot open file"
//            << exit(FatalIOError);
//    }
//    else if (!io.readHeader(isPtr()))
//    {
//        FatalIOErrorInFunction(isPtr())
//            << "problem while reading header for object " << io.name()
//            << exit(FatalIOError);
//    }
//
//    return isPtr;
//}
Foam::autoPtr<Foam::Istream> Foam::fileServers::masterFileServer::readStream
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

    fileNameList filePaths(Pstream::nProcs());
    filePaths[Pstream::myProcNo()] = fName;
    Pstream::gatherList(filePaths);

    PstreamBuffers pBufs(Pstream::nonBlocking);

    autoPtr<Istream> isPtr;

    if (Pstream::master())
    {
        //const bool uniform = uniformFile(filePaths);

        isPtr.reset(new IFstream(fName));

        // Read header
        if (!io.readHeader(isPtr()))
        {
            FatalIOErrorInFunction(isPtr())
                << "problem while reading header for object " << io.name()
                << exit(FatalIOError);
        }

DebugVar(io.headerClassName());

        word baseName;
        if
        (
            masterCollatingOFstream::isCollatingClassName
            (
                io.headerClassName(),
                baseName
            )
        )
        {
            dictionary headerDict;
            headerDict.add("version", isPtr().version());
            headerDict.add("format", isPtr().format());
            headerDict.add("class", baseName);
            if (io.note().size())
            {
                headerDict.add("note", io.note());
            }
            headerDict.add
            (
                "location",
                io.instance()/io.db().dbDir()/io.local()
            );
            headerDict.add("object", io.name());

DebugVar(headerDict);


            dictionary dict;
            dict.read(isPtr(), true);

            {
                // Extract processor0 subdict. Put in pBufs
                const word procName
                (
                    "processor" + Foam::name(Pstream::myProcNo())
                );
                const dictionary& procDict = dict.subDict(procName);

//DebugVar(procDict);


                OStringStream oss;
                // Write header
                oss.indent();
                oss.write(word("FoamFile"));
                headerDict.write(oss);
                // Write contents
                procDict.write(oss, false);

                UOPstream os(Pstream::myProcNo(), pBufs);
                string s(oss.str());
DebugVar(s);
                os.write(&s[0], s.size());
            }


            // Extract slave processor subdicts. Put in pBufs
            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
                const word procName("processor" + Foam::name(proci));
                const dictionary& procDict = dict.subDict(procName);

//DebugVar(procDict);
                OStringStream oss;
                // Write header
                oss.indent();
                oss.write(word("FoamFile"));
                headerDict.write(oss);
                // Write contents
                procDict.write(oss, false);

                UOPstream os(proci, pBufs);
                string s(oss.str());
DebugVar(s);
                os.write(&s[0], s.size());
            }

            // Close collection file (since information already extracted
            // out to pBufs)
            isPtr.clear();
        }
        else
        {
            // Already read myself. Read slave files.
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


    // isPtr will be valid on master if the file is not a Collection

    if (!isPtr.valid())
    {
        UIPstream is(Pstream::masterNo(), pBufs);
        string buf(recvSizes[Pstream::masterNo()], '\0');
        is.read(&buf[0], recvSizes[Pstream::masterNo()]);

        if (IFstream::debug)
        {
            Pout<< "Done reading " << buf.size() << " bytes" << endl;
        }

        isPtr.reset(new IStringStream(buf));

        if (!io.readHeader(isPtr()))
        {
            FatalIOErrorInFunction(isPtr())
                << "problem while reading header for object " << io.name()
                << exit(FatalIOError);
        }
    }
    return isPtr;
}


bool Foam::fileServers::masterFileServer::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    mkDir(io.path());
    fileName pathName(io.objectPath());

    autoPtr<Ostream> osPtr
    (
        new masterCollatingOFstream
        (
            io.type(),
            pathName,
            fmt,
            ver,
            cmp
        )
    );
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


Foam::autoPtr<Foam::Istream>
Foam::fileServers::masterFileServer::NewIFstream(const fileName& filePath) const
{
    if (Pstream::parRun())
    {
        // Insert logic of filePath. We assume that if a file is absolute
        // on the master it is absolute also on the slaves etc.

        fileNameList filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = filePath;
        Pstream::gatherList(filePaths);

        PstreamBuffers pBufs(Pstream::nonBlocking);

        if (Pstream::master())
        {
            const bool uniform = uniformFile(filePaths);

            if (uniform)
            {
                if (IFstream::debug)
                {
                    Pout<< "Opening global file " << filePath << endl;
                }

                // get length of file:
                off_t count(Foam::fileSize(filePath));

                std::ifstream is(filePath);
                // get length of file:
                //is.seekg(0, ios_base::end);
                //std::streamoff count = is.tellg();
                //is.seekg(0, ios_base::beg);

                if (IFstream::debug)
                {
                    Pout<< "From " << filePath
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
            string buf(recvSizes[Pstream::masterNo()], '\0');
            is.read(&buf[0], recvSizes[Pstream::masterNo()]);

            if (IFstream::debug)
            {
                Pout<< "Done reading " << buf.size() << " bytes" << endl;
            }

            // Note: IPstream is not an IStream so use a IStringStream to
            //       convert the buffer. Note that we construct with a string
            //       so it holds a copy of the buffer.
            return autoPtr<Istream>(new IStringStream(buf));
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
