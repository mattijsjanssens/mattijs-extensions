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

\*---------------------------------------------------------------------------*/

#include "commentDirective.H"
#include "dictionary.H"
#include "stringOps.H"
#include "addToMemberFunctionSelectionTable.H"
#include "dummyPrimitiveEntry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        commentDirective,
        execute,
        dictionaryIstream,
        comment
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        commentDirective,
        execute,
        primitiveEntryIstream,
        comment
    );

} // End namespace functionEntry
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::string Foam::functionEntries::commentDirective::evaluate
(
    const dictionary& parentDict,
    Istream& is
)
{
    string comment;

    auto* iss = isA<ISstream>(is);

    const_cast<ISstream&>(*iss).getLine(comment, '\n');

    return string("//") + comment;
    //return comment;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::commentDirective::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    string comment(evaluate(parentDict, is));
    //std::cout<< "Read #comment:" << comment << std::endl;
    word key("comment_" + is.name() + ':' + Foam::name(is.lineNumber()));
    //std::cout<< "key:" << key << std::endl;
    parentDict.add(new dummyPrimitiveEntry(key, comment));

    return true;
}


bool Foam::functionEntries::commentDirective::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    string comment(evaluate(parentDict, is));
//    IStringStream result(evaluate(parentDict, is));
//    entry.read(parentDict, result);
    //std::cout<< "Read #comment:" << comment << std::endl;

    // Generate unique ID
    word key("comment_" + is.name() + ':' + Foam::name(is.lineNumber()));
    //std::cout<< "key:" << key << std::endl;
    entry = dummyPrimitiveEntry(key, comment);

    return true;
}

// ************************************************************************* //
