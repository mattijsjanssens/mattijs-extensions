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
    Version of masterFileOperation that collates regIOobjects into a container
    in the processors/ subdirectory.

See also
    masterFileOperation

SourceFiles
    masterCollatingFileOperation.C

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
protected:

    // Protected classes

        //- Search for filePath; return info on how was found
        static fileName filePath(const bool, const IOobject&, pathType&, word&);


public:

        //- Runtime type information
        TypeName("masterCollatingFileOperation");


    // Constructors

        //- Construct null
        masterCollatingFileOperation();


    //- Destructor
    virtual ~masterCollatingFileOperation();


    // Member Functions

        // (reg)IOobject functionality

            //- Generate disk file name for object. Opposite of filePath.
            virtual fileName objectPath(const IOobject& io) const;

            //- Search for an object. checkGlobal : also check undecomposed case
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&
            ) const;

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
