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
    Foam::functionEntries::precompile

Description

See also

SourceFiles
    precompileEntry.C

\*---------------------------------------------------------------------------*/

#ifndef precompileEntry_H
#define precompileEntry_H

#include "ifeqEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class precompileEntry Declaration
\*---------------------------------------------------------------------------*/

class precompileEntry
:
    public ifeqEntry
{
private:

    // Private Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute
        (
            DynamicList<filePos>& stack,
            dictionary& parentDict,
            Istream& is
        );

        //- No copy construct
        precompileEntry(const precompileEntry&) = delete;

        //- No copy assignment
        void operator=(const precompileEntry&) = delete;


public:

    //- Runtime type information
    ClassName("precompile");

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
