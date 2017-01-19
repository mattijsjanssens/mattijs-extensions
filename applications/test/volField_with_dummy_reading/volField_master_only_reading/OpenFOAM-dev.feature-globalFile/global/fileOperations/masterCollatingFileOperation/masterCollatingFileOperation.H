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
    Foam::fileOperations::masterCollatingFileOperation

Description
    fileOperations that performs all file operations on the master processor.
    Requires the calls to be parallel synchronised!

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_masterCollatingFileOperation_H
#define fileOperations_masterCollatingFileOperation_H

#include "OSspecific.H"
#include "masterFileOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                Class masterCollatingFileOperation Declaration
\*---------------------------------------------------------------------------*/

class masterCollatingFileOperation
:
    public masterFileOperation
{


public:

        //- Runtime type information
        TypeName("masterCollatingFileOperation");


    // Constructors

        //- Construct null
        masterCollatingFileOperation();


    //- Destructor
    virtual ~masterCollatingFileOperation();


    // Member Functions

        // OSSpecific equivalents


        // (reg)IOobject functinality

            //- Reads header for regIOobject and returns an Istream
            //  to read the contents.
            virtual autoPtr<Istream> readStream
            (
                regIOobject&,
                const fileName&
            ) const;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state.
            virtual bool writeObject
            (
                const regIOobject&,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED
            ) const;

//            //- Generate an Istream that reads a file
//            virtual autoPtr<Istream> NewIFstream(const fileName&) const;

            //- Generate an Ostream that writes a file
            virtual autoPtr<Ostream> NewOFstream
            (
                const fileName& pathname,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
