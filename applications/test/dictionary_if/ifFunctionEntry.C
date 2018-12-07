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

#include "ifFunctionEntry.H"
#include "addToMemberFunctionSelectionTable.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "Time.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(ifFunctionEntry, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        ifFunctionEntry,
        execute,
        dictionaryIstream,
        if
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionEntries::ifFunctionEntry::readToken(token& t, Istream& is)
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


void Foam::functionEntries::ifFunctionEntry::expand
(
    const dictionary& dict,
    string& keyword
)
{
    if (keyword[0] == '$')
    {
        stringOps::inplaceExpand(keyword, dict, true, false);
    }
}


void Foam::functionEntries::ifFunctionEntry::expand
(
    const dictionary& dict,
    token& t
)
{
    if (t.isWord())
    {
        word& w = const_cast<word&>(t.wordToken());
        expand(dict, w);

        // Re-form as a string token so we can compare to string
        t = token(string(w), t.lineNumber());
    }
    else if (t.isString())
    {
        string& w = const_cast<string&>(t.stringToken());
        expand(dict, w);
    }
}


void Foam::functionEntries::ifFunctionEntry::skipUntil
(
    const dictionary& parentDict,
    const word& endWord,
    Istream& is
)
{
    while (!is.eof())
    {
        token t;
        readToken(t, is);
        if (t.isWord())
        {
            if (t.wordToken() == "#if")
            {
                skipUntil(parentDict, "#endif", is);
            }
            else if (t.wordToken() == endWord)
            {
                return;
            }
        }
    }

    FatalIOErrorInFunction(parentDict)
        << "Did not find #endif" << exit(FatalIOError);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifFunctionEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const label start = is.lineNumber();

    // Read first token and expand any string
    token cond1(is);
    expand(parentDict, cond1);

    // Read second token and expand any string
    token cond2(is);
    expand(parentDict, cond2);

    bool equal = (cond1 == cond2);

    if (equal)
    {
        Info
            << "Starting #" << typeName << " " << cond1
            << " == " << cond2
            << " at line " << start
            << " in file " <<  parentDict.name() << endl;

        while (!is.eof())
        {
            token t;
            readToken(t, is);
            if (t.isWord() && t.wordToken() == "#if")
            {
                // Evaluate
                execute(parentDict, is);
            }
            else if
            (
                t.isWord()
             && (t.wordToken() == "#else" || t.wordToken() == "#endif")
            )
            {
                if (t.wordToken() == "#else")
                {
                    // Now skip until #endif
                    skipUntil(parentDict, "#endif", is);
                }
                break;
            }
            else
            {
                is.putBack(t);
                bool ok = entry::New(parentDict, is);
                if (!ok)
                {
                    break;
                }
            }
        }
    }
    else
    {
        Info
            << "Skipping #" << typeName << " " << cond1
            << " == " << cond2
            << " at line " << start
            << " in file " <<  parentDict.name() << endl;

        // Fast-forward to #else
        token t;
        while (!is.eof())
        {
            readToken(t, is);
            if (t.isWord() && t.wordToken() == "#if")
            {
                skipUntil(parentDict, "#endif", is);
            }
            else if
            (
                t.isWord()
             && (t.wordToken() == "#else" || t.wordToken() == "#endif")
            )
            {
                break;
            }
        }

        if (t.wordToken() == "#else")
        {
            // Evaluate until #endif
            while (!is.eof())
            {
                token t;
                readToken(t, is);
                if (t.isWord() && t.wordToken() == "#if")
                {
                    // Evaluate
                    execute(parentDict, is);
                }
                else if (t.isWord() && t.wordToken() == "#endif")
                {
                    break;
                }
                else
                {
                    is.putBack(t);

                    bool ok = entry::New(parentDict, is);
                    if (!ok)
                    {
                        break;
                    }
                }
            }
        }
    }

    return true;
}


// ************************************************************************* //
