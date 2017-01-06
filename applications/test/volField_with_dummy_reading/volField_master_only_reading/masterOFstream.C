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

#include "masterOFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::Ostream> Foam::masterOFstream::open
(
    const fileName& fName
) const
{
    mkDir(fName.path());
Pout<< "Openinig file " << fName << Foam::endl;
    autoPtr<Ostream> osPtr
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

Foam::masterOFstream::masterOFstream
(
    const fileName& pathName,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OStringStream(format, version),
    pathName_(pathName),
    compression_(compression)
{
    DebugVar(pathName_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterOFstream::~masterOFstream()
{
    if (Pstream::parRun())
    {
        List<fileName> filePaths(Pstream::nProcs());
        filePaths[Pstream::myProcNo()] = pathName_;
        Pstream::gatherList(filePaths);

        if (Pstream::master())
        {
    DebugVar(filePaths);


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
    DebugVar(uniform);


            if (uniform)
            {
                const fileName& fName = object0;
                autoPtr<Ostream> osPtr(open(fName));
                osPtr() << str();
                if (!osPtr().good())
                {
                    FatalIOErrorInFunction(osPtr())
                        << "Failed writing to " << fName << exit(FatalIOError);
                }
                return;
            }
        }

    DebugVar("here");

        // Different files
        PstreamBuffers pBufs(Pstream::nonBlocking);
    DebugVar("here");

        // Send my buffer to master
        if (!Pstream::master())
        {
            UOPstream os(Pstream::masterNo(), pBufs);
            os << str();
        }
    DebugVar("here");

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

    DebugVar(recvSizes);
        if (Pstream::master())
        {
            // Write my own data
            {
                const fileName& fName = filePaths[Pstream::myProcNo()];
                autoPtr<Ostream> osPtr(open(fName));
                osPtr() << str();
                if (!osPtr().good())
                {
                    FatalIOErrorInFunction(osPtr())
                        << "Failed writing to " << fName << exit(FatalIOError);
                }
            }
    DebugVar("here");

            for (label proci = 1; proci < Pstream::nProcs(); proci++)
            {
    DebugVar("here");
                UIPstream is(proci, pBufs);
                List<char> buf(recvSizes[proci]);

Pout<< "REading from " << proci << " size:" << buf.size() << Foam::endl;
                is.read(buf.begin(), buf.size());

    DebugVar("here");
                const fileName& fName = filePaths[proci];

                autoPtr<Ostream> osPtr(open(fName));
//                osPtr().write(buf.begin(), buf.size());

    DebugVar("here");
                if (!osPtr().good())
                {
                    FatalIOErrorInFunction(osPtr())
                        << "Failed writing to " << fName
                        << exit(FatalIOError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
