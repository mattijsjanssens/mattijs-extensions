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
#include "OSspecific.H"
#include "Time.H"
#include "IFstream.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(masterFileServer, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::masterFileServer::filePath
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

        if (isFile(objectPath))
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

                if (isFile(parentObjectPath))
                {
                    searchType = fileServer::PARENTOBJECT;
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
                        searchType = fileServer::FINDINSTANCE;
                        return fName;
                    }
                }
            }
        }

        return fileName::null;
    }
}


Foam::fileName Foam::masterFileServer::objectPath
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

Foam::masterFileServer::masterFileServer()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterFileServer::~masterFileServer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::Istream>
Foam::masterFileServer::NewIFstream(IOobject& io) const
{
    if (Pstream::parRun())
    {
        // Insert logic of filePath. We assume that if a file is absolute
        // on the master it is absolute also on the slaves etc.

        List<fileName> filePaths(Pstream::nProcs());
        {
            pathType searchType = fileServer::NOTFOUND;
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
            if (searchType == fileServer::FINDINSTANCE)
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


Foam::autoPtr<Foam::Ostream> Foam::masterFileServer::NewOFstream
(
    const fileName& pathName,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return autoPtr<Ostream>(new masterOFstream(pathName, fmt, ver, cmp));
}


bool Foam::masterFileServer::mkDir(const fileName& dir) const
{
    if (Pstream::parRun())
    {
        List<fileName> filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = dir;
        Pstream::gatherList(filePaths);

        if (Pstream::master())
        {
            for (label i = 1; i < filePaths.size(); i++)
            {
                if (filePaths[i] != filePaths[0])
                {
                    Foam::mkDir(filePaths[i]);
                }
            }
            // Or scatter back result?
            return Foam::mkDir(filePaths[0]);
        }
        else
        {
            return true;
        }
    }
    else
    {
        return Foam::mkDir(dir);
    }
}


// ************************************************************************* //
