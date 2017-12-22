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
    Foam::functionObjects::parProfiling

Group
    grpUtilitiesFunctionObjects

Description

SourceFiles
    parProfiling.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_parProfiling_H
#define functionObjects_parProfiling_H

#include "functionObject.H"
#include "Tuple2.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
typedef FixedList<Tuple2<label, scalar>, 3> statData;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class parProfiling Declaration
\*---------------------------------------------------------------------------*/

class parProfiling
:
    public functionObject
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        parProfiling(const parProfiling&);

        //- Disallow default bitwise assignment
        void operator=(const parProfiling&);


public:

    //- Runtime type information
    TypeName("parProfiling");


    // Constructors

        //- Construct from Time and dictionary
        parProfiling
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~parProfiling();


    // Member Functions

        //- Read the parProfiling data
        //virtual bool read(const dictionary&);

        //- Execute file updates
        virtual bool execute();

        //- Do nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
