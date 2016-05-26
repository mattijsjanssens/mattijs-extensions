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

Application
    foamDict

Description
    Interrogates dictionaries

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Converts old scope syntax to new syntex
word scope(const fileName& entryName)
{
    if (entryName.find(':') != string::npos)
    {
        wordList entryNames(entryName.components(':'));

        word entry(entryNames[0]);
        for (label i = 1; i < entryNames.size(); i++)
        {
            entry += word('.') + entryNames[i];
        }
        return entry;
    }
    else
    {
        return entryName;
    }
}

int main(int argc, char *argv[])
{
    argList::addNote("manipulates dictionaries");

    argList::noBanner();
    argList::validArgs.append("dictionary");
    argList::addBoolOption("keywords", "report keywords");
    argList::addOption("entry", "name", "report/select the named entry");
    //argList::addOption("set", "name", "value", "sets the named entry");
    argList::addOption("value", "value", "value for the selected entry");

    #include "setRootCase.H"

    fileName dictFileName(args.rootPath()/args.caseName()/args[1]);

    autoPtr<IFstream> dictFile(new IFstream(dictFileName));

    if (dictFile().good())
    {
        bool changed = false;

        // Read but preserve headers
        dictionary dict;
        dict.read(dictFile(), true);

        if (args.optionFound("entry"))
        {
            word entryName(scope(args.option("entry")));

            const entry* entPtr = NULL;

            entPtr = dict.lookupScopedEntryPtr
            (
                entryName,
                false,
                true            // wildcards
            );

            if (entPtr)
            {
                if (args.optionFound("value"))
                {
                    Info<< "Old value:" << endl;
                    Info<< *entPtr << endl;

                    string valStr(args.option("value"));
                    valStr = string(entPtr->keyword()) + ' ' + valStr;

                    DebugVar(valStr);
                    IStringStream str(valStr);
                    autoPtr<entry> ePtr(entry::New(str));
                    DebugVar(ePtr());

//                     primitiveEntry& e =
//                         const_cast<primitiveEntry&>
//                         (
//                             dynamic_cast<const primitiveEntry&>
//                             (
//                                 *entPtr
//                             )
//                         );
// 
//                     e.read(dict, str);
// 
//                     Info<< "New value:" << endl;
//                     Info<< e << endl;
// 
//                     changed = true;

                }
                else if (args.optionFound("keywords"))
                {
                    const dictionary& dict = entPtr->dict();
                    forAllConstIter(dictionary, dict, iter)
                    {
                        Info<< iter().keyword() << endl;
                    }
                }
                else
                {
                    Info<< *entPtr << endl;
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot find entry "
                    << entryName
                    << " in dictionary " << dictFileName;
                FatalError.exit(2);
            }
        }
        else if (args.optionFound("keywords"))
        {
            forAllConstIter(dictionary, dict, iter)
            {
                Info<< iter().keyword() << endl;
            }
        }
        else
        {
            Info<< dict;
        }

        if (changed)
        {
            dictFile.clear();
            OFstream os(dictFileName);
            dict.write(os, false);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot open file " << dictFileName;
        FatalError.exit(1);
    }

    return 0;
}


// ************************************************************************* //
