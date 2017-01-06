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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "IFstream.H"

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
autoPtr<Istream> createStream(const IOobject& io)
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
Pout<< "Master detected type:" << label(searchType)
    << " at:" << filePaths[Pstream::myProcNo()] << endl;

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
Pout<< "Slave constructed filename:" << filePaths[Pstream::myProcNo()] << endl;
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
autoPtr<Istream> readStream(IOobject& io)
{
    autoPtr<Istream> isPtr(createStream(io));

    if (!io.readHeader(isPtr()))
    {
        FatalIOErrorInFunction(isPtr())
            << "problem while reading header for object " << io.name()
            << exit(FatalIOError);
    }

    return isPtr;
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

    {
        autoPtr<Istream> isPtr(readStream(io));
        IOdictionary dict(io, isPtr());
        DebugVar(dict);
    }

    return 0;
}


// ************************************************************************* //
