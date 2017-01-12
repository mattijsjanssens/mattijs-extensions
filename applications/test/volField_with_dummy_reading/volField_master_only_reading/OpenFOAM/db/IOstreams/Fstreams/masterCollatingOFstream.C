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

#include "masterCollatingOFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamBuffers.H"
#include "IOdictionary.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::OFstream> Foam::masterCollatingOFstream::open
(
    const fileName& fName
) const
{
    mkDir(fName.path());
    autoPtr<OFstream> osPtr
    (
        new OFstream
        (
            fName,
            IOstream::BINARY,   //format(),
            version(),
            compression_
        )
    );
    if (!osPtr().good())
    {
        FatalIOErrorInFunction(osPtr())
            << "Could not open file " << fName
            << exit(FatalIOError);
    }
    return osPtr;
}


std::streamoff Foam::masterCollatingOFstream::writeBuffers
(
    List<std::streamoff>& start,
    List<std::streamoff>& size
) const
{
    start.setSize(Pstream::nProcs());
    size.setSize(Pstream::nProcs());

    if (Pstream::parRun())
    {
        List<fileName> filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = pathName_;
        Pstream::gatherList(filePaths);

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
                const fileName& fName = object0;

                autoPtr<OFstream> osPtr(open(fName));
                osPtr().writeQuoted(str(), false);
                if (!osPtr().good())
                {
                    FatalIOErrorInFunction(osPtr())
                        << "Failed writing to " << fName << exit(FatalIOError);
                }
                return osPtr().stdStream().tellp();
            }
        }

        // Different files
        PstreamBuffers pBufs(Pstream::nonBlocking);

        // Send my buffer to master
        if (!Pstream::master())
        {
            UOPstream os(Pstream::masterNo(), pBufs);
            string s(this->str());
            os.write(&s[0], s.size());
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (Pstream::master())
        {
            const fileName& fName = filePaths[Pstream::myProcNo()];
            autoPtr<OFstream> osPtr(open(fName));

            // We don't have IOobject so cannot use writeHeader
            OFstream& os = osPtr();
            IOobject::writeBanner(os)
                << "FoamFile\n{\n"
                << "    version     " << version() << ";\n"
                << "    format      " << format() << ";\n"
                << "    class       " << collatingClassName(typeName_) << ";\n"
                << "    location    " << fName << ";\n"
                << "    object      " << fName.name() << ";\n"
                << "}" << nl;
            IOobject::writeDivider(os) << nl;


            // Write my own data
            {
                const word procName
                (
                    "processor"
                  + Foam::name(Pstream::myProcNo())
                );
                os  << procName << " {\n";

                start[Pstream::myProcNo()] = os.stdStream().tellp();

                os.writeQuoted(str(), false);
                if (!os.good())
                {
                    FatalIOErrorInFunction(os)
                        << "Failed writing to " << fName << exit(FatalIOError);
                }

                size[Pstream::myProcNo()] =
                    os.stdStream().tellp()
                   -start[Pstream::myProcNo()];

                os  << "}\n\n";
            }

            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
                UIPstream is(proci, pBufs);
                List<char> buf(recvSizes[proci]);

                is.read(buf.begin(), buf.size());

                const fileName& fName = filePaths[proci];

                const word procName("processor" + Foam::name(proci));
                os  << procName<< " {\n";

                start[proci] = os.stdStream().tellp();

                os.writeQuoted(string(buf.begin(), buf.size()), false);
                if (!os.good())
                {
                    FatalIOErrorInFunction(os)
                        << "Failed writing to " << fName
                        << exit(FatalIOError);
                }

                size[proci] = os.stdStream().tellp()-start[proci];

                os  << "}\n\n";
            }
            IOobject::writeEndDivider(os);

            return os.stdStream().tellp();
        }
        else
        {
            return 0;
        }
    }
    else
    {
        autoPtr<OFstream> osPtr(open(pathName_));

        start[Pstream::myProcNo()] = 0;

        osPtr().writeQuoted(str(), false);

        size[Pstream::myProcNo()] =
            osPtr().stdStream().tellp()
           -start[Pstream::myProcNo()];

        if (!osPtr().good())
        {
            FatalIOErrorInFunction(osPtr())
                << "Failed writing to " << pathName_ << exit(FatalIOError);
        }

        return osPtr().stdStream().tellp();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterCollatingOFstream::masterCollatingOFstream
(
    const fileName& typeName,
    const fileName& pathName,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OStringStream(format, version),
    typeName_(typeName),
    pathName_(pathName),
    compression_(compression)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterCollatingOFstream::~masterCollatingOFstream()
{
    List<std::streamoff> start;
    List<std::streamoff> size;
    std::streamoff fileSize = writeBuffers(start, size);

    dictionary dict;
    dict.add("size", fileSize);

    forAll(start, proci)
    {
        Info<< "    processor " << proci << " start:" << start[proci]
            << " size:" << size[proci] << Foam::endl;

        dictionary procDict;
        // TBD: 64 bits!
        procDict.add("start", start[proci]);
        procDict.add("size", size[proci]);
        dict.add(word("processor" + Foam::name(proci)), procDict);
    }

    // Write as a dictionary
    OFstream os(pathName_ + ".index");
    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Problem writing index file" << exit(FatalIOError);
    }

    IOobject::writeBanner(os)
        << "FoamFile\n{\n"
        << "    version     " << version() << ";\n"
        << "    format      " << format() << ";\n"
        << "    class       " << IOdictionary::typeName << ";\n"
        << "    location    " << os.name() << ";\n"
        << "    object      " << os.name().name() << ";\n"
        << "}" << nl;
    IOobject::writeDivider(os) << nl;

    dict.write(os, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::masterCollatingOFstream::collatingClassName
(
    const word& baseName
)
{
    return word("Collection<" + baseName + ">");
}


bool Foam::masterCollatingOFstream::isCollatingClassName
(
    const word& className,
    word& baseName
)
{
    size_t off = className.find("Collection<");

    if (off != string::npos)
    {
        baseName = className.substr(11, className.size()-11);

        label sz = baseName.size();
        if (baseName[sz-1] == '>')
        {
            baseName.resize(sz-1);
            return true;
        }
    }
    return false;
}


std::streamoff Foam::masterCollatingOFstream::readIndexFile
(
    const fileName& pathName,
    List<std::streamoff>& start,
    List<std::streamoff>& size
)
{
    std::streamoff overallSize = 0;

    IFstream is(pathName + ".index");
    if (!is.good())
    {
        start.clear();
        size.clear();
    }
    else
    {
        dictionary dict(is);

        dict.readIfPresent("size", overallSize);

        start.setSize(Pstream::nProcs(), -1);
        size.setSize(Pstream::nProcs(), -1);

        forAll(start, proci)
        {
            const word procName("processor" + Foam::name(proci));

            if (dict.found(procName))
            {
                const dictionary& procDict = dict.subDict(procName);

                procDict.readIfPresent("start", start[proci]);
                procDict.readIfPresent("size", size[proci]);
            }
        }
    }

    return overallSize;
}


// ************************************************************************* //
