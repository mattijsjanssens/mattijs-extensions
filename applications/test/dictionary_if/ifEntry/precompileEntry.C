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

#include "IOstream.H"
#include "precompileEntry.H"
#include "addToMemberFunctionSelectionTable.H"
#include "IStringStream.H"
#include "Switch.H"
#include "IOstreams.H"
#include <iostream>
#include "int64.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(precompileEntry, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        precompileEntry,
        execute,
        dictionaryIstream,
        precompile
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//bool Foam::functionEntries::precompileEntry::execute
//(
//    DynamicList<filePos>& stack,
//    dictionary& parentDict,
//    Istream& is
//)
//{
//    const label nNested = stack.size();
//
//    stack.append(filePos(is.name(), is.lineNumber()));
//
//    // Read line
//    string line;
//    dynamic_cast<ISstream&>(is).getLine(line);
//    line += ';';
//    IStringStream lineStream(line);
//    const primitiveEntry e("precompileEntry", parentDict, lineStream);
//    const Switch doIf(e.stream());
//
//    Info
//        << "Using #" << typeName << " " << doIf
//        << " at line " << stack.last().second()
//        << " in file " <<  stack.last().first() << endl;
//
//    bool ok = ifeqEntry::execute(doIf, stack, parentDict, is);
//
//    if (stack.size() != nNested)
//    {
//        FatalIOErrorInFunction(parentDict)
//            << "Did not find matching #endif for condition starting"
//            << " at line " << stack.last().second()
//            << " in file " <<  stack.last().first() << exit(FatalIOError);
//    }
//
//    return ok;
//}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::precompileEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
DebugVar(is.info());
DebugVar(parentDict);
//DebugVar(is.type());
    if (isA<ISstream>(is))
    {
        ISstream& iss = dynamic_cast<ISstream&>(is);
        std::istream& stdis = iss.stdStream();

        std::istream::streampos here = stdis.tellg();
        DebugVar(here);

        DebugVar(stdis.good());
        DebugVar(stdis.fail());
        DebugVar(stdis.bad());


        dictionary save(parentDict);

//        executedictionaryIstreamMemberFunctionTable dictTable
//        (
//            *executedictionaryIstreamMemberFunctionTablePtr_
//        );
//        DebugVar(dictTable.sortedToc());


        executeprimitiveEntryIstreamMemberFunctionTable& et =
            *executeprimitiveEntryIstreamMemberFunctionTablePtr_;

        // Preserve old entries
        executeprimitiveEntryIstreamMemberFunctionTable entryTable(et);

        // Replace entries
        forAllConstIter
        (
            preProcesspreProcessMemberFunctionTable,
            *preProcesspreProcessMemberFunctionTablePtr_,
            iter
        )
        {
            et.set(iter.key(), iter());
        }


        // Read dictionary
        parentDict.read(is);

        // Restore stream

        //- Use rewind
        //is.rewind();

        //- Use direct positioning
        //stdis.seekg(here, stdis.beg);
        //stdis.setstate(iostate::good)

        //- Use direct positioning
        //stdis.rdbuf()->pubseekpos(here);
        Pout<< "Before seek at:" << int64_t(stdis.tellg()) << endl;
        stdis.seekg(here, stdis.beg);
        Pout<< "Now again at:" << int64_t(stdis.tellg()) << endl;


        //stdis.clear(stdis.eofbit);
        //stdis.clear(stdis.failbit);
        //stdis.clear(stdis.badbit);
        stdis.clear(std::ios::failbit);
        //stdis.setstate(stdis.rdstate());

        DebugVar(stdis.good());
        DebugVar(stdis.fail());
        DebugVar(stdis.bad());

        // Take over std state to ISstream level
        is.setState(stdis.rdstate());

        //stdis.clear();
        //is.setGood();

        DebugVar(is.good());
        DebugVar(is.fail());
        DebugVar(is.bad());


DebugVar(is.info());
        //is.print(Pout);

        // Restore table
        et.transfer(entryTable);

        // Restore dictionary
        parentDict = save;

        // Re-read
DebugVar(is.info());
        //stdis.rdbuf()->pubseekpos(here);
        //is.setState(stdis.rdstate());
Pout<< nl << "** STARTING REREADING" << endl;
        dictionary d2(is);
Pout<< nl << "** DONE REREADING d2:" << d2 << endl;
DebugVar(is.info());
    }
    return true;
}


// ************************************************************************* //
