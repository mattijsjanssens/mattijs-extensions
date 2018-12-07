/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::functionEntries::ifEq

Description

SourceFiles
    ifEq.C

\*---------------------------------------------------------------------------*/

#ifndef ifEq_H
#define ifEq_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class dlLibraryTable;

namespace functionEntries
{

// Forward declaration of friend classes
class calcEntry;

/*---------------------------------------------------------------------------*\
                         Class ifEq Declaration
\*---------------------------------------------------------------------------*/

class ifEq
:
    public functionEntry
{
    // Private Member Functions

        static void readToken(token& t, Istream& is);

        //- No copy construct
        ifEq(const ifEq&) = delete;

        //- No copy assignment
        void operator=(const ifEq&) = delete;


public:

    //- Runtime type information
    ClassName("ifEq");


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);

//         //- Execute the functionEntry in a primitiveEntry context
//         static bool execute
//         (
//             const dictionary& parentDict,
//             primitiveEntry& entry,
//             Istream& is
//         );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
