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

#ifndef ifFunctionEntry_H
#define ifFunctionEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                       Class ifFunctionEntry Declaration
\*---------------------------------------------------------------------------*/

class ifFunctionEntry
:
    public functionEntry
{
    // Private Member Functions

        //- Read tokens. Skip dummy tokens
        static void readToken(token& t, Istream& is);

        //- In-place expand a string
        static void expand
        (
            const dictionary& dict,
            string& keyword
        );

        //- In-place expand a string or word token
        static void expand
        (
            const dictionary& dict,
            token& t
        );

        //- No copy construct
        ifFunctionEntry(const ifFunctionEntry&) = delete;

        //- No copy assignment
        void operator=(const ifFunctionEntry&) = delete;


public:

    //- Runtime type information
    ClassName("if");


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
