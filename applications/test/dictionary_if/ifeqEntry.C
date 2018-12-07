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

#include "ifeqEntry.H"
#include "addToMemberFunctionSelectionTable.H"
#include "IStringStream.H"
#include "stringOps.H"
#include "ifEntry.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(ifeqEntry, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        ifeqEntry,
        execute,
        dictionaryIstream,
        ifeq
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionEntries::ifeqEntry::readToken(token& t, Istream& is)
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


void Foam::functionEntries::ifeqEntry::expand
(
    const dictionary& dict,
    string& keyword,
    token& t
)
{
    if (keyword[0] == '$')
    {
        word varName = keyword(1, keyword.size()-1);

        // lookup the variable name in the given dictionary
        const entry* ePtr = dict.lookupScopedEntryPtr
        (
            varName,
            true,
            true
        );
        if (ePtr)
        {
            t = token(ePtr->stream());
        }
        else
        {
            // String expansion
            stringOps::inplaceExpand(keyword, dict, true, false);
            // Re-form as a string token so we can compare to string
            t = token(keyword, t.lineNumber());
        }
    }
    else if (!t.isString())
    {
        // Re-form as a string token so we can compare to string
        t = token(keyword, t.lineNumber());
    }
}


void Foam::functionEntries::ifeqEntry::expand
(
    const dictionary& dict,
    token& t
)
{
    if (t.isWord())
    {
        expand(dict, const_cast<word&>(t.wordToken()), t);
    }
    else if (t.isVariable())
    {
        expand(dict, const_cast<string&>(t.stringToken()), t);
    }
    else if (t.isString())
    {
        expand(dict, const_cast<string&>(t.stringToken()), t);
    }
}


void Foam::functionEntries::ifeqEntry::skipUntil
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
            if (t.wordToken() == "#if" || t.wordToken() == "#ifeq")
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
        << "Did not find matching " << endWord << exit(FatalIOError);
}


bool Foam::functionEntries::ifeqEntry::execute
(
    const bool doIf,
    const label start,
    dictionary& parentDict,
    Istream& is
)
{
    if (doIf)
    {
        while (!is.eof())
        {
            token t;
            readToken(t, is);

            if (t.isWord() && t.wordToken() == "#ifeq")
            {
                // Recurse to evaluate
                execute(parentDict, is);
            }
            else if (t.isWord() && t.wordToken() == "#if")
            {
                // Recurse to evaluate
                ifEntry::execute(parentDict, is);
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
        // Fast-forward to #else
        token t;
        while (!is.eof())
        {
            readToken(t, is);
            if
            (
                t.isWord()
             && (t.wordToken() == "#if" || t.wordToken() == "#ifeq")
            )
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
                if (t.isWord() && t.wordToken() == "#ifeq")
                {
                    // Recurse to evaluate
                    execute(parentDict, is);
                }
                else if (t.isWord() && t.wordToken() == "#if")
                {
                    // Recurse to evaluate
                    ifEntry::execute(parentDict, is);
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifeqEntry::execute
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

    const bool equal = (cond1 == cond2);

    Info
        << "Evaluating #" << typeName << " " << cond1
        << " == " << cond2
        << " at line " << start
        << " in file " <<  parentDict.name() << endl;

    return execute(equal, start, parentDict, is);
}


// ************************************************************************* //
