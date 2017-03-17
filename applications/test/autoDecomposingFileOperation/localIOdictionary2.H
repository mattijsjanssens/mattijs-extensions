/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::localIOdictionary2

Description
    localIOdictionary2 is derived from IOdictionary but excludes parallel
    master reading.

SourceFiles
    localIOdictionary2.C

\*---------------------------------------------------------------------------*/

#ifndef localIOdictionary2_H
#define localIOdictionary2_H

#include "baseIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localIOdictionary2 Declaration
\*---------------------------------------------------------------------------*/

class localIOdictionary2
:
    public baseIOdictionary
{

public:

    // Constructors

        //- Construct given an IOobject
        localIOdictionary2(const IOobject& io);

        //- Construct given an IOobject and dictionary
        localIOdictionary2(const IOobject& io, const dictionary& dict);

        //- Construct given an IOobject and Istream
        localIOdictionary2(const IOobject& io, Istream& is);

        //- Construct given an IOobject and wanted typeName
        localIOdictionary2(const IOobject& io, const word& wantedType);


    //- Destructor
    virtual ~localIOdictionary2();


    // Member functions

        //- Is object global
        virtual bool global() const
        {
            return false;
        }

        //- Return complete path + object name if the file exists
        //  in the case otherwise null
        virtual fileName filePath() const
        {
            // Use default (local only) search strategy
            return localFilePath();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
