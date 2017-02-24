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
#include "OFstream.H"
#include "OSspecific.H"
#include <pthread.h>
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void* Foam::threadedOFstream::writeFile(void *threadarg)
{
    //Pout<< "threadedOFstream::writeFile" << Foam::endl;
    writeData* dataPtr = static_cast<writeData*>(threadarg);
    writeData& data = *dataPtr;

    //Pout<< "threadedOFstream::writeFile:"
    //    << " writing " << data.data_.size()
    //    << " bytes to file " << data.pathName_ << Foam::endl;

    mkDir(data.pathName_.path());
    OFstream os
    (
        data.pathName_,
        data.format_,
        data.version_,
        data.compression_
    );

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Could not open file"
            << exit(FatalIOError);
    }
    os.writeQuoted(data.data_, false);

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Failed writing data of size " << data.data_.size()
            << exit(FatalIOError);
    }

    //Pout<< "threadedOFstream::writeFile:"
    //    << " succesfully wrote " << data.pathName_ << Foam::endl;

    delete dataPtr;

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threadedOFstream::threadedOFstream
(
    pthread_t& thread,
    const fileName& pathName,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OStringStream(format, version),
    thread_(thread),
    pathName_(pathName),
    compression_(compression)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::threadedOFstream::~threadedOFstream()
{
    writeData* dataPtr
    (
        new writeData
        (
            pathName_,
            compression(),
            format(),
            version(),
            str()
        )
    );

    int rc = pthread_create
    (
        &thread_,
        nullptr,
        writeFile,
        dataPtr
    );

    if (rc)
    {
        FatalErrorInFunction << "problem creating thread" << exit(FatalError);
    }
}


// ************************************************************************* //
