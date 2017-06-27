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
    Foam::IStringStream2

Description
    Input from memory buffer stream.

SourceFiles
    StringStreamsPrint.C

\*---------------------------------------------------------------------------*/

#ifndef IStringStream2_H
#define IStringStream2_H

#include "ISstream.H"
#include <sstream>
#include "zstr.hpp"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class IStringStream2 Declaration
\*---------------------------------------------------------------------------*/

class IStringStream2
:
    public ISstream
{
    // Private data

        //- Data
        std::stringbuf buf_;


public:

    // Constructors

        //- Construct from string
        IStringStream2
        (
            const string& str,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType cmp=UNCOMPRESSED
        )
        :
            ISstream
            (
                (
                    cmp == compressionType::UNCOMPRESSED
                 ? *reinterpret_cast<std::istream*>
                    (
                        new std::istream(&buf_)
                    )
                 : *reinterpret_cast<std::istream*>
                    (
                        new zstr::istream(&buf_)
                    )
                ),
                "IStringStream2.sourceFile",
                format,
                version,
                cmp
            ),
            buf_(str)
        {}


        //- Construct from char*
        IStringStream2
        (
            const char* buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType cmp=UNCOMPRESSED
        )
        :
            ISstream
            (
                (
                    cmp == compressionType::UNCOMPRESSED
                  ? *reinterpret_cast<std::istream*>(new std::istream(&buf_))
                  : *reinterpret_cast<std::istream*>(new zstr::istream(&buf_))
                ),
                "IStringStream2.sourceFile",
                format,
                version
            ),
            buf_(buffer)
        {}


    //- Destructor
    ~IStringStream2()
    {
        if (compression() == compressionType::UNCOMPRESSED)
        {
            delete &dynamic_cast<std::istream&>(stdStream());
        }
        else
        {
            delete &dynamic_cast<zstr::istream&>(stdStream());
        }
    }


    // Member functions

        // Access

            //- Return the string
            string str() const
            {
                return buf_.str();
            }


        // Print

            //- Print description to Ostream
            void print(Ostream&) const;


    // Member operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary:
        //  e.g. thing thisThing(IFstream("thingFileName")());
        Istream& operator()() const
        {
            return const_cast<IStringStream2&>(*this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
