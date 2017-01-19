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
#include "IOobject.H"
#include "decomposedBlockData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterCollatingOFstream::masterCollatingOFstream
(
    const fileName& pathName,
    const UPstream::commsTypes commsType,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OStringStream(format, version),
    pathName_(pathName),
    commsType_(commsType),
    compression_(compression)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterCollatingOFstream::~masterCollatingOFstream()
{
    autoPtr<OSstream> osPtr;
    if (UPstream::master())
    {
        // Note: always write binary. These are strings so readable
        //       anyway. They have already be tokenised on the sending side.
        osPtr.reset
        (
            new OFstream
            (
                pathName_,
                IOstream::BINARY,
                version(),
                compression()
            )
        );

        //writeHeader(osPtr());
        // We don't have IOobject so cannot use writeHeader
        OSstream& os = osPtr();
        IOobject::writeBanner(os)
            << "FoamFile\n{\n"
            << "    version     " << version() << ";\n"
            << "    format      " << os.format() << ";\n"
            << "    class       " << decomposedBlockData::typeName << ";\n"
            << "    location    " << pathName_ << ";\n"
            << "    object      " << pathName_.name() << ";\n"
            << "}" << nl;
        IOobject::writeDivider(os) << nl;
    }


    string s(str());
    UList<char> slice(const_cast<char*>(s.data()), label(s.size()));

    List<std::streamoff> start;
    decomposedBlockData::writeBlocks(osPtr, start, slice, commsType_);
    if (osPtr.valid() && !osPtr().good())
    {
        FatalIOErrorInFunction(osPtr())
            << "Failed writing to " << pathName_ << exit(FatalIOError);
    }
}


// ************************************************************************* //
