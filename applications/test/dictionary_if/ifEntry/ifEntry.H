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
    Foam::functionEntries::if

Description
    Conditional parsing of dictionary entries.

    E.g.
    \verbatim
        U_inlet 15;

        #if #calc "${U_inlet} < 10"
            ..
        #else
            ..
        #endif
    \endverbatim

    Note:
    - only supports single line - no '\' is supported
    - condition should be readable as a Switch (supports 0,1, true, false, etc.)

See also
    Foam::functionEntries::ifeq
    Foam::Switch

SourceFiles
    ifEntry.C

\*---------------------------------------------------------------------------*/

#ifndef ifEntry_H
#define ifEntry_H

#include "ifeqEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class ifEntry Declaration
\*---------------------------------------------------------------------------*/

class ifEntry
:
    public ifeqEntry
{
private:

    friend class ifeqEntry;

    // Private Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute
        (
            DynamicList<filePos>& stack,
            dictionary& parentDict,
            Istream& is
        );

        //- No copy construct
        ifEntry(const ifEntry&) = delete;

        //- No copy assignment
        void operator=(const ifEntry&) = delete;


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
