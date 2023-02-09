/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Application
    dictionaryTest

Description

\*---------------------------------------------------------------------------*/

#include "ISstream.H"
#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "dictionary.H"
#include "stringOps.H"
#include "Function1.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void expand(const dictionary& topDict, dictionary& d)
{
    for (auto& e : d)
    {
        if (e.isStream())
        {
            const tokenList& toks = e.stream();
            for (const token& t : toks)
            {
                if (t.isWord())
                {
                    string& s = const_cast<string&>(t.stringToken());
                    Pout<< "    string:" << s << " size:" << s.size()
                        << endl;
                    if (s.size() > 2 && s[0] == '\\' && s[1] == '$')
                    {
                        // Remove leading '\\'
                        s = s.substr(1);

                        Pout<< "Expanding var:" << s << endl;

                        stringOps::inplaceExpand
                        (
                            s,
                            topDict,
                            true,
                            false,
                            false
                        );
                    }
                }
            }
        }
        else
        {
            const dictionary& subDict = e.dict();
            expand(topDict, const_cast<dictionary&>(subDict));
        }
    }
}


//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addArgument("dict .. dictN");
    argList args(argc, argv, false, true);

//    {
//        dictionary dict;
//        dict.add(word("ab" + getEnv("WM_MPLIB") + "cd"), 16);
//        Info<< "dict:" << dict << endl;
//
//        string s("DDD_${ab${WM_MPLIB}cd}_EEE");
//        stringOps::inplaceExpand(s, dict, true, false);
//        Info<< "variable expansion:" << s << endl;
//    }

    // 1. Determine syntax for unparsed variable name
    //dictionary dict;
    //dict.add(new primitiveEntry("name", "\\$WM_MPLIB"));

    Pout<<" ** STARTING **" << endl;
    ISstream::keepComments = true;
    dictionary dict(IFstream("testDict")());
Pout << dict << endl;
    Pout<<" ** STOPPING **" << endl;
    ISstream::keepComments = false;

return 0;


    //- not possible - needs Field<word> possible
    //auto fnPtr = Function1<word>::New("bla", dict);


    // 2. One-step replacement of variable name
    expand(dict, dict);
    DebugVar(dict);

    // 3. Or could do two step:
    //  - rewrite string containing '\\$' into 'variable' token
    //  - expand variable token



    //const auto* ePtr = dict.findEntry("name");
    //auto& is = ePtr->stream();
    //DebugVar(is.size());
    //DebugVar(is.name());

//    token& t = is[0];
//    t = token(string("$WM_MPLIB"));
//    t.setType(token::VARIABLE);
//    DebugVar(t.isVariable());
//    DebugVar(t);

    return 0;
}


// ************************************************************************* //
