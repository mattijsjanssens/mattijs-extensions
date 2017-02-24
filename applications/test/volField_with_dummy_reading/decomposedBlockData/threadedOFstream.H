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

Class
    Foam::threadedOFstream

Description
    Threaded drop-in replacement for OFstream.

SourceFiles
    threadedOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef threadedOFstream_H
#define threadedOFstream_H

#include "OStringStream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class threadedOFstream Declaration
\*---------------------------------------------------------------------------*/

class threadedOFstream
:
    public OStringStream
{
    // Private class

        class writeData
        {
        public:

            fileName pathName_;
            IOstream::compressionType compression_;
            IOstream::streamFormat format_;
            IOstream::versionNumber version_;
            string data_;

            writeData
            (
                const fileName& pathName,
                IOstream::compressionType compression,
                IOstream::streamFormat format,
                IOstream::versionNumber version,
                const string& data
            )
            :
                pathName_(pathName),
                compression_(compression),
                format_(format),
                version_(version),
                data_(data)
            {}
        };


    // Private data

        pthread_t& thread_;

        const fileName pathName_;

        const IOstream::compressionType compression_;


    // Private Member Functions

        static void* writeFile(void *threadarg);


public:

    // Constructors

        //- Construct and set stream status. Thread is only initialised
        //  at destruction.
        threadedOFstream
        (
            pthread_t& thread,
            const fileName& pathname,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    //- Destructor
    ~threadedOFstream();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
