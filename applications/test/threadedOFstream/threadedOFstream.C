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

#include "threadedOFstream.H"
#include "OFstreamWriter.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threadedOFstream::threadedOFstream
(
    OFstreamWriter& writer,
    const fileName& pathName,
    streamFormat format,
    versionNumber version,
    compressionType compression,
    const bool append
)
:
    OStringStream(format, version),
    writer_(writer),
    pathName_(pathName),
    compression_(compression),
    append_(append)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::threadedOFstream::~threadedOFstream()
{
    bool ok = writer_.write
    (
        pathName_,
        str(),
        format(),
        version(),
        compression_,   // use supplied, not the one from OStringStream
        append_
    );

    if (!ok)
    {
        FatalErrorInFunction
            << "Did not write file " << pathName_ << exit(FatalError);
    }
}


// ************************************************************************* //
