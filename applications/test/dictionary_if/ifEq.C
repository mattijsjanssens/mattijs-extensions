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

#include "ifEq.H"
#include "addToMemberFunctionSelectionTable.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(ifEq, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        ifEq,
        execute,
        dictionaryIstream,
        ifEq
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionEntries::ifEq::readToken(token& t, Istream& is)
{
    // Skip dummy tokens - avoids entry::getKeyword consuming #else, #endif
    do
    {
        if
        (
            is.read(t).bad()
         || is.eof()
         || !t.good()
        )
        {
            return;
        }
    }
    while (t == token::END_STATEMENT);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifEq::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const token cond1(is);
    const token cond2(is);

    if (cond1 == cond2)
    {
        const label start = is.lineNumber();

        Info
            << "Starting #" << typeName << " at line " << start
            << " in file " <<  parentDict.name() << endl;

        while (!is.eof())
        {
            token t;
            readToken(t, is);
            if (t.isWord() && t.wordToken() == "#else")
            {
                const label line = is.lineNumber();

                // Now skip until #endif
                while (!is.eof())
                {
                    token t;
                    readToken(t, is);
                    if (t.isWord() && t.wordToken() == "#endif")
                    {
                        Info
                            << "Skipped #else from line " << line
                            << " to " << is.lineNumber()
                            << " in file " <<  parentDict.name() << endl;
                        break;
                    }
                }
                break;
            }
            is.putBack(t);

            bool ok = entry::New(parentDict, is);
            if (!ok)
            {
                break;
            }
        }
    }
    else
    {
        // Fast-forward to #else
        while (!is.eof())
        {
            token t;
            readToken(t, is);
            if (t.isWord() && t.wordToken() == "#else")
            {
                break;
            }
        }

        const label line = is.lineNumber();
        Info
            << "Starting #else at line " << line
            << " in file " <<  parentDict.name() << endl;

        // Parse until #endif
        while (!is.eof())
        {
            token t;
            readToken(t, is);
            if (t.isWord() && t.wordToken() == "#endif")
            {
                Info
                    << "Finished #else from line " << line
                    << " to " << is.lineNumber()
                    << " in file " <<  parentDict.name() << endl;
                break;
            }
            is.putBack(t);

            bool ok = entry::New(parentDict, is);
            if (!ok)
            {
                break;
            }
        }
    }

    return true;
}


// ************************************************************************* //
