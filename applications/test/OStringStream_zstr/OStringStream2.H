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
    Foam::OStringStream2

Description
    Output to memory buffer stream.

SourceFiles
    StringStreamsPrint.C

\*---------------------------------------------------------------------------*/

#ifndef OStringStream2_H
#define OStringStream2_H

#include "OSstream.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class OStringStream2 Declaration
\*---------------------------------------------------------------------------*/

class OStringStream2
:
    public OSstream
{

public:

    // Public classes

        //- Less function class that can be used for sorting
        class zstr_ostringstream
        :
            public zstr::ostream
        {
        public:

            std::stringbuf buf_;

            zstr_ostringstream()
            :
                zstr::ostream(&buf_)
            {}

            string str() const
            {
                return buf_.str();
            }
        };


    // Constructors

        //- Construct and set stream status
        OStringStream2
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType cmp=UNCOMPRESSED
        )
        :
            OSstream
            (
                (
                    cmp == compressionType::UNCOMPRESSED
                  ? *(new std::ostringstream())
                  : *(new zstr_ostringstream())
                ),
                "OStringStream.sinkFile",
                format,
                version,
                cmp
            )
        {}

//        //- Construct as copy
//        OStringStream2(const OStringStream2& oss)
//        :
//            OSstream
//            (
//               *(
//                    new std::ostringstream
//                    (
//                        dynamic_cast<const std::ostringstream&>
//                        (
//                            oss.stdStream()
//                        ).str()
//                    )
//                ),
//                oss.name(),
//                oss.format(),
//                oss.version()
//            )
//        {}


    //- Destructor
    ~OStringStream2()
    {
        if (compression() == compressionType::UNCOMPRESSED)
        {
            delete &dynamic_cast<std::ostringstream&>(stdStream());
        }
        else
        {
            delete &dynamic_cast<zstr_ostringstream&>(stdStream());
        }
    }


    // Member functions

        // Access

            //- Return the string
            string str() const
            {
                if (compression() == compressionType::UNCOMPRESSED)
                {
                    return dynamic_cast<const std::ostringstream&>
                    (
                        stdStream()
                    ).str();
                }
                else
                {
                    return dynamic_cast<const zstr_ostringstream&>
                    (
                        stdStream()
                    ).str();
                }
            }


        // Edit

            //- Rewind the OStringStream
            void rewind()
            {
                stdStream().rdbuf()->pubseekpos(0);
            }


        // Print

            //- Print description to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
