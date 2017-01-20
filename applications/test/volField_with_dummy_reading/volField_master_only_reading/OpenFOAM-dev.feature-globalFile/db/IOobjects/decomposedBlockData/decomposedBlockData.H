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
    Foam::decomposedBlockData

Description
    decomposedBlockData is a List<char> with IO on the master processor only.

SourceFiles
    decomposedBlockData.C

\*---------------------------------------------------------------------------*/

#ifndef decomposedBlockData_H
#define decomposedBlockData_H

#include "IOList.H"
#include "regIOobject.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class decomposedBlockData Declaration
\*---------------------------------------------------------------------------*/

class decomposedBlockData
:
    public regIOobject,
    public List<char>
{
protected:

    // Protected data

        //- Type to use for gather
        const UPstream::commsTypes commsType_;


public:

    TypeName("decomposedBlockData");


    // Constructors

        //- Construct given an IOobject
        decomposedBlockData
        (
            const IOobject&,
            const UPstream::commsTypes = UPstream::scheduled
        );

        //- Construct given an IOobject and List<char>
        decomposedBlockData
        (
            const IOobject&,
            const List<char>&,
            const UPstream::commsTypes = UPstream::scheduled
        );

        //- Construct by transferring the IOList contents
        decomposedBlockData
        (
            const IOobject&,
            const Xfer<List<char>>&,
            const UPstream::commsTypes = UPstream::scheduled
        );


    //- Destructor
    virtual ~decomposedBlockData();


    // Member functions

        //- Read object
        virtual bool read();

        //- Pure virtual writaData function.
        //  Must be defined in derived types
        virtual bool writeData(Ostream&) const
        {
            NotImplemented;
            return false;
        }

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType
        ) const;


    // Helpers

        //- Read data into *this. Istream is only valid on master.
        static bool readBlocks
        (
            autoPtr<Istream>& isPtr,
            List<char>& data,
            const UPstream::commsTypes commsType
        );

        //- Write *this. Ostream only valid on master. Returns starts of
        //  processor blocks
        static bool writeBlocks
        (
            autoPtr<OSstream>& osPtr,
            List<std::streamoff>& start,
            const UList<char>&,
            const UPstream::commsTypes
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
