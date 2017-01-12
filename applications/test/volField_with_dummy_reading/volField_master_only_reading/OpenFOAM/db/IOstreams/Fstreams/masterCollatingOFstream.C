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
                return;
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
                  + Foam::name(Pstream::masterNo())
                );
                os  << procName << " {\n";

                Info<< "Starting entry for processor " << procName
                    << " at " << os.stdStream().tellp() << Foam::endl;

                os.writeQuoted(str(), false);
                if (!os.good())
                {
                    FatalIOErrorInFunction(os)
                        << "Failed writing to " << fName << exit(FatalIOError);
                }

                Info<< "Finished entry processor " << procName
                    << " at " << os.stdStream().tellp()
                    << Foam::endl;

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

                Info<< "Starting entry for processor " << procName
                    << " at " << os.stdStream().tellp() << Foam::endl;

                os.writeQuoted(string(buf.begin(), buf.size()), false);
                if (!os.good())
                {
                    FatalIOErrorInFunction(os)
                        << "Failed writing to " << fName
                        << exit(FatalIOError);
                }

                Info<< "Finished entry processor " << procName
                    << " at " << os.stdStream().tellp()
                    << Foam::endl;

                os  << "}\n\n";
            }
            IOobject::writeEndDivider(os);
        }
    }
    else
    {
        autoPtr<OFstream> osPtr(open(pathName_));
        osPtr().writeQuoted(str(), false);
        if (!osPtr().good())
        {
            FatalIOErrorInFunction(osPtr())
                << "Failed writing to " << pathName_ << exit(FatalIOError);
        }
    }
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


// ************************************************************************* //
