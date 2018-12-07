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
    Foam::functionEntries::ifeq

Description
    Conditional parsing of dictionary entries.

    E.g.
    \verbatim
        a #calc "0.123";
        b 1.23e-1;

        #ifeq $a $b
            ..
        #else
            ..
        #endif
    \endverbatim

See also
    Foam::functionEntries::if

SourceFiles
    ifeqEntry.C

\*---------------------------------------------------------------------------*/

#ifndef ifeqEntry_H
#define ifeqEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class ifeqEntry Declaration
\*---------------------------------------------------------------------------*/

class ifeqEntry
:
    public functionEntry
{
protected:

    // Protected Member Functions

        //- Read tokens. Skip dummy tokens
        static void readToken(token& t, Istream& is);

        //- In-place expand a variable (string/word/var starting with '$')
        static void expand
        (
            const dictionary& dict,
            string& keyword,
            token& t
        );

        //- In-place expand a string or word token
        static void expand
        (
            const dictionary& dict,
            token& t
        );

        //- Consume tokens until reached a specific word
        static void skipUntil
        (
            const dictionary& parentDict,
            const word& endWord,
            Istream& is
        );

        //- Main driver: depending on 'equal' starts evaluating or
        //  skips forward to #else
        static bool execute
        (
            const bool equal,
            const label start,
            dictionary& parentDict,
            Istream& is
        );

private:

    // Private Member Functions

        //- No copy construct
        ifeqEntry(const ifeqEntry&) = delete;

        //- No copy assignment
        void operator=(const ifeqEntry&) = delete;


public:

    //- Runtime type information
    ClassName("ifeq");

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
