/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::functionEntries::commentDirective

Description
    Place holder to insert '//' comment lines. Behaves like #message but
    prints with '//'.

SourceFiles
    commentDirective.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntries_commentDirective_H
#define functionEntries_commentDirective_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                      Class commentDirective Declaration
\*---------------------------------------------------------------------------*/

class commentDirective
:
    public functionEntry
{
    // Private Member Functions

    //- Evaluate
    static string evaluate(const dictionary& parentDict, Istream& is);


public:

    //- Execute in a primitiveEntry context.
    //  Reports message string (after expansion) - does not alter entry.
    static bool execute
    (
        const dictionary& parentDict,
        primitiveEntry& entry,
        Istream& is
    );

    //- Execute in a sub-dict context.
    //  Reports message string (after expansion) - does not alter dictionary.
    static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
