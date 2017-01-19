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

#include "masterFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "IFstream.H"
#include "masterOFstream.H"
#include "masterCollatingOFstream.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(masterFileOperation, 0);
    addToRunTimeSelectionTable(fileOperation, masterFileOperation, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::masterFileOperation::filePath
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
        fileName path = io.path();
        fileName objectPath = path/io.name();

        if (Foam::isFile(objectPath))
        {
            searchType = fileOperation::OBJECT;
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
                    searchType = fileOperation::PARENTOBJECT;
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
                        searchType = fileOperation::FINDINSTANCE;
                        return fName;
                    }
                }
            }
        }

        return fileName::null;
    }
}


Foam::fileName Foam::fileOperations::masterFileOperation::objectPath
(
    const IOobject& io,
    const pathType& searchType,
    const word& instancePath
)
{
    // Replacement for IOobject::objectPath()

    switch (searchType)
    {
        case fileOperation::ABSOLUTE:
        {
            return io.instance()/io.name();
        }
        break;

        case fileOperation::OBJECT:
        {
            return io.path()/io.name();
        }
        break;

        case fileOperation::PARENTOBJECT:
        {
            return
                io.rootPath()/io.time().globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case fileOperation::FINDINSTANCE:
        {
            return
                io.rootPath()/io.caseName()
               /instancePath/io.db().dbDir()/io.local()/io.name();
        }
        break;

        case fileOperation::NOTFOUND:
        {
            return fileName::null;
        }
        break;
    }
}


bool Foam::fileOperations::masterFileOperation::uniformFile
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


bool Foam::fileOperations::masterFileOperation::collating() const
{
    if (collating_ == -1)
    {
        bool collate = debug::optimisationSwitch("collateFiles", false);
        if (collate)
        {
            WarningInFunction << "Switching on parallel file writing collation"
                << endl;
            collating_ = 1;
        }
        else
        {
            WarningInFunction << "No parallel file writing collation" << endl;
            collating_ = 0;
        }
    }

    return (collating_ == 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::masterFileOperation::masterFileOperation()
:
    collating_(-1)
{
    if (regIOobject::fileModificationChecking == regIOobject::timeStampMaster)
    {
        WarningInFunction
            << "Resetting fileModificationChecking to timeStamp" << endl;
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
    }
    else if
    (
        regIOobject::fileModificationChecking
     == regIOobject::inotifyMaster
    )
    {
        WarningInFunction
            << "Resetting fileModificationChecking to inotifyMaster" << endl;
        regIOobject::fileModificationChecking = regIOobject::inotifyMaster;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::masterFileOperation::~masterFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperations::masterFileOperation::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return masterOp<mode_t, mkDirOp>(dir, mkDirOp(mode));
}


bool Foam::fileOperations::masterFileOperation::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return masterOp<mode_t, chModOp>(fName, chModOp(mode));
}


mode_t Foam::fileOperations::masterFileOperation::mode
(
    const fileName& fName
) const
{
    return masterOp<mode_t, modeOp>(fName, modeOp());
}


Foam::fileName::Type Foam::fileOperations::masterFileOperation::type
(
    const fileName& fName
) const
{
    return fileName::Type(masterOp<label, typeOp>(fName, typeOp()));
}


bool Foam::fileOperations::masterFileOperation::exists
(
    const fileName& fName,
    const bool checkGzip
) const
{
    return masterOp<bool, existsOp>(fName, existsOp(checkGzip));
}


bool Foam::fileOperations::masterFileOperation::isDir
(
    const fileName& fName
) const
{
    return masterOp<bool, isDirOp>(fName, isDirOp());
}


bool Foam::fileOperations::masterFileOperation::isFile
(
    const fileName& fName,
    const bool checkGzip
) const
{
    return masterOp<bool, isFileOp>(fName, isFileOp());
}


off_t Foam::fileOperations::masterFileOperation::fileSize
(
    const fileName& fName
) const
{
    return masterOp<off_t, fileSizeOp>(fName, fileSizeOp());
}


time_t Foam::fileOperations::masterFileOperation::lastModified
(
    const fileName& fName
) const
{
    return masterOp<time_t, lastModifiedOp>(fName, lastModifiedOp());
}


double Foam::fileOperations::masterFileOperation::highResLastModified
(
    const fileName& fName
) const
{
    return masterOp<double, lastModifiedHROp>
    (
        fName,
        lastModifiedHROp()
    );
}


bool Foam::fileOperations::masterFileOperation::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return masterOp<bool, mvBakOp>(fName, mvBakOp(ext));
}


bool Foam::fileOperations::masterFileOperation::rm(const fileName& fName) const
{
    return masterOp<bool, rmOp>(fName, rmOp());
}


bool Foam::fileOperations::masterFileOperation::rmDir
(
    const fileName& dir
) const
{
    return masterOp<bool, rmDirOp>(dir, rmDirOp());
}


Foam::fileNameList Foam::fileOperations::masterFileOperation::readDir
(
    const fileName& dir,
    const fileName::Type type,
    const bool filtergz
) const
{
    return masterOp<fileNameList, readDirOp>
    (
        dir,
        readDirOp(type, filtergz)
    );
}


bool Foam::fileOperations::masterFileOperation::cp
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterOp<bool, cpOp>(src, dst, cpOp());
}


bool Foam::fileOperations::masterFileOperation::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterOp<bool, lnOp>(src, dst, lnOp());
}


bool Foam::fileOperations::masterFileOperation::mv
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterOp<bool, mvOp>(src, dst, mvOp());
}


Foam::fileName Foam::fileOperations::masterFileOperation::filePath
(
    const IOobject& io
) const
{
    fileName objPath;
    pathType searchType = fileOperation::NOTFOUND;
    word newInstancePath;
    if (Pstream::master())
    {
        objPath = filePath(io, searchType, newInstancePath);
    }
    label masterType(searchType);
    Pstream::scatter(masterType);
    searchType = pathType(masterType);
    if (searchType == fileOperation::FINDINSTANCE)
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
Foam::fileOperations::masterFileOperation::objectStream
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


//Foam::autoPtr<Foam::Istream>
//Foam::fileOperations::masterFileOperation::readStream
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
//            "masterFileOperation::readStream()",
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
Foam::autoPtr<Foam::Istream>
Foam::fileOperations::masterFileOperation::readStream
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

        autoPtr<IFstream> ifsPtr(new IFstream(fName));
        IFstream& is = ifsPtr();

        // Read header
        if (!io.readHeader(is))
        {
            FatalIOErrorInFunction(is)
                << "problem while reading header for object " << io.name()
                << exit(FatalIOError);
        }

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
            // Construct dictionary with header information so we can
            // stream this to the individual processors
            dictionary headerDict;
            {
                headerDict.add("version", is.version());
                headerDict.add("format", is.format());
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
                //DebugVar(headerDict);
            }


            // See if index file present
            List<std::streamoff> start;
            List<std::streamoff> size;
            std::streamoff overallSize = masterCollatingOFstream::readIndex
            (
                fName,
                start,
                size
            );

            forAll(size, proci)
            {
                if (size[proci] > labelMax)
                {
                    FatalIOErrorInFunction(fName)
                        << "Size of slice for processor " << proci
                        << " exceeds range of label " << labelMax
                        << exit(FatalIOError);
                }
            }

            const std::streamoff realSize = Foam::fileSize(fName);

            if (overallSize && overallSize != realSize)
            {
                IOWarningInFunction(is)
                    << "Ignoring index file since contains wrong size "
                    << overallSize << " and file size is " << realSize
                    << endl;
                overallSize = 0;
            }


            if (overallSize && overallSize == realSize)
            {
                // Have index file. Can slice file directly
                if (IFstream::debug)
                {
                    Pout<< "For file " << fName
                        << " found index file with size:" << overallSize
                        << " and starts:" << start << Foam::endl;
                }

                string headerString;
                {
                    // Get header as string
                    OStringStream oss;
                    oss.indent();
                    oss.write(word("FoamFile"));
                    headerDict.write(oss);
                    headerString = oss.str();
                }


                {
                    // Read master slice from file
                    List<char> buf
                    (
                        static_cast<label>(size[Pstream::masterNo()])
                    );
                    is.stdStream().seekg
                    (
                        start[Pstream::masterNo()],
                        ios_base::beg
                    );
                    is.stdStream().read
                    (
                        buf.begin(),
                        size[Pstream::masterNo()]
                    );

                    // Prepend header
                    string s(headerString);
                    label sz = s.size();
                    s.resize(sz+buf.size());
                    forAll(buf, i)
                    {
                        s[sz++] = buf[i];
                    }
//DebugVar(s);
                    UOPstream os(Pstream::myProcNo(), pBufs);
                    os.write(&s[0], s.size());
                }
                for (label proci = 1; proci < Pstream::nProcs(); proci++)
                {
                    // Read slice from file

                    List<char> buf(static_cast<label>(size[proci]));
                    is.stdStream().seekg(start[proci], ios_base::beg);
                    is.stdStream().read(buf.begin(), size[proci]);

                    // Prepend header
                    string s(headerString);
                    label sz = s.size();
                    s.resize(sz+buf.size());
                    forAll(buf, i)
                    {
                        s[sz++] = buf[i];
                    }
//DebugVar(s);
                    UOPstream os(proci, pBufs);
                    os.write(&s[0], s.size());
                }
            }
            else
            {
                // Read complete dictionary
                dictionary dict;
                dict.read(is, true);

                OStringStream oss;
                {
                    // Extract processor0 subdict. Put in pBufs
                    const word procName
                    (
                        "processor" + Foam::name(Pstream::myProcNo())
                    );
                    const dictionary& procDict = dict.subDict(procName);

                    // Write header
                    oss.indent();
                    oss.write(word("FoamFile"));
                    headerDict.write(oss);
                    // Write contents
                    procDict.write(oss, false);
                    string s(oss.str());
//DebugVar(s);
                    UOPstream os(Pstream::myProcNo(), pBufs);
                    os.write(&s[0], s.size());
                }


                // Extract slave processor subdicts. Put in pBufs
                for (label proci = 1; proci < Pstream::nProcs(); proci++)
                {
                    const word procName("processor" + Foam::name(proci));
                    const dictionary& procDict = dict.subDict(procName);

                    oss.rewind();
                    // Write header
                    oss.indent();
                    oss.write(word("FoamFile"));
                    headerDict.write(oss);
                    // Write contents
                    procDict.write(oss, false);
                    string s(oss.str());
//DebugVar(s);
                    UOPstream os(proci, pBufs);
                    os.write(&s[0], s.size());
                }
            }

            // Close collection file (since information already extracted
            // out to pBufs)
            ifsPtr.clear();
        }
        else
        {
            // Open master (steal from ifsPtr)
            isPtr.reset(ifsPtr.ptr());

            // Read slave files
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
                List<char> buf(static_cast<label>(count));
                is.read(buf.begin(), count);

                UOPstream os(proci, pBufs);
                os.write(buf.begin(), count);
            }
        }
    }

    labelList recvSizes;
    pBufs.finishedSends(recvSizes);

    // isPtr will be valid on master if the file is not a Collection. In
    // all cases the information is in the PstreamBuffers

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


bool Foam::fileOperations::masterFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    mkDir(io.path());
    fileName pathName(io.objectPath());

    autoPtr<Ostream> osPtr;

    //if (collating())
    //{
    //    osPtr.reset
    //    (
    //        new masterCollatingOFstream
    //        (
    //            io.type(),
    //            pathName,
    //            fmt,
    //            ver,
    //            cmp
    //        )
    //    );
    //}
    //else
    //{
        osPtr.reset
        (
            new masterOFstream
            (
                pathName,
                fmt,
                ver,
                cmp
            )
        );
    //}
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
Foam::fileOperations::masterFileOperation::NewIFstream
(
    const fileName& filePath
) const
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

                if (IFstream::debug)
                {
                    Pout<< "From " << filePath
                        <<  " reading " << label(count) << " bytes" << endl;
                }
                List<char> buf(static_cast<label>(count));
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

                    off_t count(Foam::fileSize(filePaths[proci]));

                    std::ifstream is(filePaths[proci]);

                    if (IFstream::debug)
                    {
                        Pout<< "From " << filePaths[proci]
                            <<  " reading " << label(count) << " bytes" << endl;
                    }
                    List<char> buf(static_cast<label>(count));
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


Foam::autoPtr<Foam::Ostream>
Foam::fileOperations::masterFileOperation::NewOFstream
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
