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
    Foam::masterCollatingOFstream

Description
    Master-only drop-in replacement for OFstream.

SourceFiles
    masterCollatingOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef masterCollatingOFstream_H
#define masterCollatingOFstream_H

#include "OStringStream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstream;

/*---------------------------------------------------------------------------*\
                       Class masterCollatingOFstream Declaration
\*---------------------------------------------------------------------------*/

class masterCollatingOFstream
:
    public OStringStream
{
    // Private data

        const fileName typeName_;

        const fileName pathName_;

        const IOstream::compressionType compression_;


    // Private Member Functions

        //- Open file with checking
        autoPtr<OFstream> open(const fileName&) const;


public:

    // Constructors

        //- Construct and set stream status
        masterCollatingOFstream
        (
            const fileName& headerClassName,
            const fileName& pathname,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    //- Destructor
    ~masterCollatingOFstream();


    // Member Functions

        //- Construct headerClassName given typeName
        static word collatingClassName(const word&);

        //- Is className a collating type? Extract typeName if so
        static bool isCollatingClassName(const word&, word&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
