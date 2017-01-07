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
    Test-volField

Description
    The idea is that argList allocates a fileServer (through command-line
    arguments or environment vars). This has operations to do
    - file existence checking
    - mkDir, rm etc.
    - open an IFstream
    - open an OFstream
    and everywhere instead of constructing an I/OFstream we do a
        fileServer::NewIFstream(io)
        fileServer::NewOFstream(io)

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "IFstream.H"
#include "OFstream.H"
#include "masterOFstream.H"
#include "fileServer.H"
#include "masterFileServer.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

enum pathType
{
    NOTFOUND,           // not found
    ABSOLUTE,           // instance is absolute directory
    OBJECT,             // objectPath exists
    PARENTOBJECT,       // parent of object path
    FINDINSTANCE        // file found in time directory
};


// Replacement for IOobject::filePath()
fileName filePath
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
        if (isFile(objectPath))
        {
            searchType = ABSOLUTE;
            return objectPath;
        }
        else
        {
            searchType = NOTFOUND;
            return fileName::null;
        }
    }
    else
    {
        fileName path = io.path();
        fileName objectPath = path/io.name();

        if (isFile(objectPath))
        {
            searchType = OBJECT;
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

                if (isFile(parentObjectPath))
                {
                    searchType = PARENTOBJECT;
                    return parentObjectPath;
                }
            }

            if (!isDir(path))
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

                    if (isFile(fName))
                    {
                        searchType = FINDINSTANCE;
                        return fName;
                    }
                }
            }
        }

        return fileName::null;
    }
}


// Replacement for IOobject::objectPath()
fileName objectPath
(
    const IOobject& io,
    const pathType& searchType,
    const word& instancePath
)
{
    switch (searchType)
    {
        case pathType::ABSOLUTE:
        {
            return io.instance()/io.name();
        }
        break;

        case pathType::OBJECT:
        {
            return io.path()/io.name();
        }
        break;

        case pathType::PARENTOBJECT:
        {
            return
                io.rootPath()/io.time().globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case pathType::FINDINSTANCE:
        {
            return
                io.rootPath()/io.caseName()
               /instancePath/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case pathType::NOTFOUND:
        {
            return fileName::null;
        }
        break;
    }
}


// Master-only reading of files
autoPtr<Istream> createReadStream(const IOobject& io)
{
    if (Pstream::parRun())
    {
        // Insert logic of filePath. We assume that if a file is absolute
        // on the master it is absolute also on the slaves etc.

        List<fileName> filePaths(Pstream::nProcs());
        {
            pathType searchType = NOTFOUND;
            word newInstancePath;
            if (Pstream::master())
            {
                filePaths[Pstream::myProcNo()] = filePath
                (
                    io,
                    searchType,
                    newInstancePath
                );
            }
            label masterType(searchType);
            Pstream::scatter(masterType);
            searchType = pathType(masterType);
            if (searchType == pathType::FINDINSTANCE)
            {
                Pstream::scatter(newInstancePath);
            }

            if (!Pstream::master())
            {
                filePaths[Pstream::myProcNo()] = objectPath
                (
                    io,
                    searchType,
                    newInstancePath
                );
            }
        }
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
                Pout<< "Reading " << io.objectPath()
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
        pathType searchType;
        word newInstancePath;
        fileName fName(filePath(io, searchType, newInstancePath));

        return autoPtr<Istream>(new IFstream(fName));
    }
}


// Replacement for regIOobject::readStream
autoPtr<Istream> readStream(const fileServer& server, IOobject& io)
{
    if (IFstream::debug)
    {
        InfoInFunction
            << "Reading object " << io.name()
            << " from file " << io.objectPath()
            << endl;
    }

    if (io.readOpt() == IOobject::NO_READ)
    {
        FatalErrorInFunction
            << "NO_READ specified for read-constructor of object " << io.name()
            << " of class " << io.headerClassName()
            << abort(FatalError);
    }

    // Construct object stream and read header if not already constructed
    //if (!isPtr_)
    autoPtr<Istream> isPtr;
    {

        fileName objPath;
        //if (io.watchIndex() != -1)
        //{
        //    // File is being watched. Read exact file that is being watched.
        //    objPath = io.time().getFile(io.watchIndex());
        //}
        //else
        {
            // Search intelligently for file
            objPath = server.filePath(io);

            if (!objPath.size())
            {
                FatalIOError
                (
                    "regIOobject::readStream()",
                    __FILE__,
                    __LINE__,
                    io.objectPath(),
                    0
                )   << "cannot find file"
                    << exit(FatalIOError);
            }
        }

        //if (!(isPtr_ = objectStream(objPath)))
        isPtr = server.NewIFstream(objPath);
        if (!io.readHeader(isPtr()))
        {
            FatalIOErrorInFunction(isPtr())
                << "problem while reading header for object " << io.name()
                << exit(FatalIOError);
        }
    }

    //// Mark as uptodate if read successfully
    //if (io.watchIndex() != -1)
    //{
    //    io.time().setUnmodified(io.watchIndex());
    //}

    return isPtr;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Master-only writing
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Replacement for regIOobject::writeObject
bool writeObject
(
    const fileServer& server,
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
)
{
    if (!io.good())
    {
        SeriousErrorInFunction
            << "bad object " << io.name()
            << endl;

        return false;
    }

    if (io.instance().empty())
    {
        SeriousErrorInFunction
            << "instance undefined for object " << io.name()
            << endl;

        return false;
    }

    if
    (
        io.instance() != io.time().timeName()
     && io.instance() != io.time().system()
     && io.instance() != io.time().caseSystem()
     && io.instance() != io.time().constant()
     && io.instance() != io.time().caseConstant()
    )
    {
        const_cast<regIOobject&>(io).instance() = io.time().timeName();
    }

    if (OFstream::debug)
    {
        InfoInFunction << "Writing file " << io.objectPath();
    }


    bool osGood = false;

    {
        // Try opening an OFstream for object
        //mkDir(io.path());
        //OFstream os(io.objectPath(), fmt, ver, cmp);
        autoPtr<Ostream> osPtr
        (
            server.NewOFstream
            (
                io.objectPath(),
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

        io.writeEndDivider(os);

        osGood = os.good();
    }

    if (OFstream::debug)
    {
        Info<< " .... written" << endl;
    }

    // Only update the lastModified_ time if this object is re-readable,
    // i.e. lastModified_ is already set
    if (io.watchIndex() != -1)
    {
        io.time().setUnmodified(io.watchIndex());
    }

    return osGood;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

//     IOobject io
//     (
//         "p",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     );
    IOobject io
    (
        "fvSolution",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    //masterFileServer server;
    autoPtr<fileServer> serverPtr
    (
        fileServer::New
        (
            fileServers::masterFileServer::typeName
        )
    );
    const fileServer& server = serverPtr();

    {
        autoPtr<Istream> isPtr(readStream(server, io));
        IOdictionary dict(io, isPtr());
        DebugVar(dict);

        runTime++;
        dict.instance() = runTime.timeName();
        writeObject
        (
            server,
            dict,
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED
        );
    }

    return 0;
}


// ************************************************************************* //
