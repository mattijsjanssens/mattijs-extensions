/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::OStringStream

Description
    Output to memory buffer stream.

SourceFiles
    StringStreamsPrint.C

\*---------------------------------------------------------------------------*/

#ifndef OStringStream_H
#define OStringStream_H

#include "OSstream.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class OStringStream Declaration
\*---------------------------------------------------------------------------*/

class OStringStream
:
    public OSstream
{

    static std::ostringstream* allocator(const label size)
    {
        string s;
        s.reserve(size);

        std::cout<< "Allocating from size:" << size << std::endl;

        return new std::ostringstream(s);
    }

public:

    // Constructors

        //- Construct and set stream status
        OStringStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OSstream
            (
               *(new std::ostringstream()),
                "OStringStream.sinkFile",
                format,
                version
            )
        {}

        //- Construct with initial size and set stream status
        OStringStream
        (
            const label size,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OSstream
            (
                *allocator(size),
                "OStringStream.sinkFile",
                format,
                version
            )
        {}

        //- Construct as copy
        OStringStream(const OStringStream& oss)
        :
            OSstream
            (
               *(
                    new std::ostringstream
                    (
                        dynamic_cast<const std::ostringstream&>
                        (
                            oss.stdStream()
                        ).str()
                    )
                ),
                oss.name(),
                oss.format(),
                oss.version()
            )
        {}


    //- Destructor
    ~OStringStream()
    {
        delete &dynamic_cast<std::ostringstream&>(stdStream());
    }


    // Member functions

        // Access

            //- Return the string
            string str() const
            {
                return dynamic_cast<const std::ostringstream&>
                (
                    stdStream()
                ).str();
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
