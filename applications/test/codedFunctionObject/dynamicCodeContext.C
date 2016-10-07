/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "OSHA1stream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCodeContext::dynamicCodeContext(const dictionary& dict)
{
    read(dict);
}


Foam::dynamicCodeContext::dynamicCodeContext()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::dynamicCodeContext::clear()
{
    filterVars_.clear();
    stream_.rewind();
    sha1_.clear();
}


// const Foam::string& Foam::dynamicCodeContext::filterVar
// (
//     const word& var,
//     const word& defaultVal
// ) const
// {
//     HashTable<string>::const_iterator fnd = filterVars_.find(var);
//     if (fnd != filterVars_.end())
//     {
//         return var;
//     }
//     else
//     {
//         return defaultVal;
//     }
// }


void Foam::dynamicCodeContext::read(const dictionary& dict)
{
    forAllConstIter(dictionary, dict, iter)
    {
        const entry& e = iter();

        if (!e.isDict())
        {
            const word& key = e.keyword();

            bool isCode = true;
            if (key == "codeOptions" || key == "codeLibs")
            {
                isCode = false;
            }

            addFilterVariable(isCode, dict, key);
        }
    }
}


void Foam::dynamicCodeContext::write(Ostream& os, const dictionary& dict) const
{
    wordList keys(filterVars_.sortedToc());

    forAll(keys, i)
    {
        const word& key = keys[i];

        if (dict.found(key))
        {
            os.writeKeyword(key)
                << token::HASH << token::BEGIN_BLOCK;

            os.writeQuoted(string(dict[key]), false)
                << token::HASH << token::END_BLOCK
                << token::END_STATEMENT << nl;
        }
    }
}


void Foam::dynamicCodeContext::addFilterVariable
(
    const bool addDirective,
    const dictionary& dict,
    const word& key
)
{
    const entry* ePtr = dict.lookupEntryPtr(key, false, false);
    if (ePtr)
    {
        string s(stringOps::trim(ePtr->stream()));
        stringOps::inplaceExpand(s, dict);

        // Add to sha1
        stream_ << s;

        // Add line directive
        
        if (addDirective)
        {
            dynamicCodeContext::addLineDirective
            (
                s,
                ePtr->startLineNumber(),
                dict.name()
            );
        }

        filterVars_.set(key, s);
    }
}


void Foam::dynamicCodeContext::addFilterVariables
(
    const fileName& srcFile,
    const dictionary& dict
)
{
    IFstream is(srcFile);

    if (!is.good())
    {
        FatalIOErrorInFunction(dict)
            << "Failed opening " << srcFile
            << exit(FatalIOError);
    }

    // Extract all variables from file
    string line;
    wordHashSet vars;
    do
    {
        is.getLine(line);
        stringOps::variables(vars, line);
    }
    while (is.good());

    forAllConstIter(wordHashSet, vars, iter)
    {
        addFilterVariable(true, dict, iter.key());
    } 
}


void Foam::dynamicCodeContext::addLineDirective
(
    string& code,
    const label lineNum,
    const fileName& name
)
{
    code = "#line " + Foam::name(lineNum + 1) + " \"" + name + "\"\n" + code;
}


// ************************************************************************* //
