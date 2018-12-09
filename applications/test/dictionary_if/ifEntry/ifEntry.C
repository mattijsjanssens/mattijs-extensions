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

\*---------------------------------------------------------------------------*/

#include "ifEntry.H"
#include "addToMemberFunctionSelectionTable.H"
#include "IStringStream.H"
#include "Switch.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(ifEntry, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        ifEntry,
        execute,
        dictionaryIstream,
        if
    );
}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    DynamicList<filePos> stack(10);
    stack.append(filePos(parentDict.name(), is.lineNumber()));

    // Read line
    string line;
    dynamic_cast<ISstream&>(is).getLine(line);
    line += ';';
    IStringStream lineStream(line);
    const primitiveEntry e("ifEntry", parentDict, lineStream);
    const Switch doIf(e.stream());

    Info
        << "Using #" << typeName << " " << doIf
        << " at line " << stack.last().second()
        << " in file " <<  stack.last().first() << endl;

    bool ok = ifeqEntry::execute(doIf, stack, parentDict, is);

DebugVar(stack);

    if (stack.size() != 1)
    {
        FatalIOErrorInFunction(parentDict)
            << "Did not find matching #endif for condition starting"
            << " at line " << stack.last().second()
            << " in file " <<  stack.last().first() << exit(FatalIOError);
    }

    return ok;
}


// ************************************************************************* //
