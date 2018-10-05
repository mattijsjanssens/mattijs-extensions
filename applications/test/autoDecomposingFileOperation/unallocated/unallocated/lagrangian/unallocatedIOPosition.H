/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::unallocatedIOPosition

Description
    Helper IO class to read and write particle positions

SourceFiles
    unallocatedIOPosition.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedIOPosition_H
#define unallocatedIOPosition_H

#include "regIOobject.H"
#include "basicParticle.H"
#include "passiveParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class unallocatedIOPosition Declaration
\*---------------------------------------------------------------------------*/

class unallocatedIOPosition
:
    public regIOobject,
    public IDLList<basicParticle>
{

    // Private data


public:

    // Static data

        //- Runtime type name information. Use cloud type.
        virtual const word& type() const
        {
            return passiveParticleCloud::typeName;
        }


    // Constructors

        //- Construct read
        unallocatedIOPosition(const IOobject& io);


    // Member functions

        virtual bool readData(Istream&);

        virtual bool write(const bool valid = true) const;

        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
