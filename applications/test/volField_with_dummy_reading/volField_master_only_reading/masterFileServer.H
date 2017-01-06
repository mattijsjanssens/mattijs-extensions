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
    Foam::masterFileServer

\*---------------------------------------------------------------------------*/

#ifndef masterFileServer_H
#define masterFileServer_H

#include "fileServer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class masterFileServer Declaration
\*---------------------------------------------------------------------------*/

class masterFileServer
:
    public fileServer
{
//    // Private data
//
//    // Private Member Functions

        //- Search for filePath; return info on how was found
        static fileName filePath(const IOobject&, pathType&, word&);

        //- Construct filePath
        static fileName objectPath
        (
            const IOobject&,
            const pathType&,
            const word&
        );


public:

        //- Runtime type information
        TypeName("masterFileServer");


    // Constructors

        //- Construct null
        masterFileServer();


    //- Destructor
    virtual ~masterFileServer();


    // Member Functions

        //- Generate an Istream that reads a file
        virtual autoPtr<Istream> NewIFstream(IOobject& io) const;

        //- Generate an Ostream that reads a file
        virtual autoPtr<Ostream> NewOFstream
        (
            const fileName& pathname,
            IOstream::streamFormat format=IOstream::ASCII,
            IOstream::versionNumber version=IOstream::currentVersion,
            IOstream::compressionType compression=IOstream::UNCOMPRESSED
        ) const;

        //- Make directory
        virtual bool mkDir(const fileName&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
