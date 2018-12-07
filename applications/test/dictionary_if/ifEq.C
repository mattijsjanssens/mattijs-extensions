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

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        ifEq,
        execute,
        primitiveEntryIstream,
        ifEq
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::ifEq::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    Info
        << "Using #ifEq at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    const word bla1(is);
    const word bla2(is);

//    if (bla1 == bla2)
//    {
//        while
//        (
//            !is.eof()
//         && entry::New(*this, is)
//        )
//        {}
//
//
//
//    // get code dictionary
//    // must reference parent for stringOps::expand to work nicely
//    dictionary codeDict("#ifEq", parentDict, is);
//
////    streamingFunctionType function = getFunction(parentDict, codeDict);
////
////    // use function to write stream
////    OStringStream os(is.format());
////    (*function)(os, parentDict);
////
////    // get the entry from this stream
////    IStringStream resultStream(os.str());
////    entry.read(parentDict, resultStream);

    return true;
}


bool Foam::functionEntries::ifEq::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    Info
        << "Using #ifEq at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    const word bla1(is);
DebugVar(bla1);
    const word bla2(is);
DebugVar(bla2);

    //const dictionary currentState(parentDict);

    if (bla1 == bla2)
    {
        while (!is.eof())
        {
            autoPtr<entry> ePtr = entry::New(is);
DebugVar(is.info());
            if (ePtr.valid())
            {
                if (ePtr().isStream())
                {
                    token t(ePtr().stream());
Pout<< "Have stream:" << t << endl;
                    if (t.isWord() && t.wordToken() == "#else")
                    {
                        return true;
                    }
                }
                else
                {
Pout<< "Have dictionary:" << ePtr().dict() << endl;
                }
                parentDict.add(ePtr());
            }
        }
    }

    return true;
}


// ************************************************************************* //
