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
    Foam::OFstreamWriter

Description
    Threaded drop-in replacement for OFstream.

SourceFiles
    OFstreamWriter.C

\*---------------------------------------------------------------------------*/

#ifndef OFstreamWriter_H
#define OFstreamWriter_H

#include "IOstream.H"
#include "List.H"
#include "FIFOStack.H"
//#include <pthread.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class OFstreamWriter Declaration
\*---------------------------------------------------------------------------*/

class OFstreamWriter
{
    // Private class

        class writeData
        {
        public:

            const fileName pathName_;
            const string data_;
            const IOstream::streamFormat format_;
            const IOstream::versionNumber version_;
            const IOstream::compressionType compression_;
            const bool append_;

            writeData
            (
                const fileName& pathName,
                const string& data,
                IOstream::streamFormat format,
                IOstream::versionNumber version,
                IOstream::compressionType compression,
                const bool append
            )
            :
                pathName_(pathName),
                data_(data),
                format_(format),
                version_(version),
                compression_(compression),
                append_(append)
            {}
        };


    // Private data

        const off_t maxBufferSize_;

        //pthread_mutex_t mutex_;
        label mutex_;

        //pthread_t thread_;
        label thread_;

        FIFOStack<writeData*> objects_;

        bool threadRunning_;


    // Private Member Functions

        //- Write actual file
        static bool writeFile
        (
            const fileName& fName,
            const string& data,
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool append
        );

        //- Write all files in stack
        static void* writeAll(void *threadarg);


public:

    // Declare name of the class and its debug switch
    TypeName("OFstreamWriter");


    // Constructors

        //- Construct from buffer size. 0 = do not use thread
        OFstreamWriter(const off_t maxBufferSize);


    //- Destructor
    virtual ~OFstreamWriter();


    // Member functions

        //- Write file with contents. Blocks until writethread has space
        //  available (total file sizes < maxBufferSize)
        bool write
        (
            const fileName&,
            const string& data,
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool append
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
