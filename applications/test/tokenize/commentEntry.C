/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "commentEntry.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(commentEntry, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        commentEntry,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        commentEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::commentEntry::commentEntry
(
    const word& key,
    const dictionary& dict,
    Istream& is
)
:
    functionEntry(key, dict, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::commentEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    word comment(is);
    DebugVar(comment);
    return true;
}


bool Foam::functionEntries::commentEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    word comment(is);
    DebugVar(comment);
    return true;
}


void Foam::functionEntries::commentEntry::write(Ostream& os) const
{
    DebugVar("here");
    // Contents should be single string token
    const token& t = operator[](0);
    const string& s = t.stringToken();

    for (size_t i = 0; i < s.size(); i++)
    {
        os.write(s[i]);
    }

    os << endl;
}


// ************************************************************************* //
