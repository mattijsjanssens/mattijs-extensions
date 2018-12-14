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
#include "OSspecific.H"

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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::precompileEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    if (isA<ISstream>(is))
    {
        Info<< "#" << typeName
            << " : pre-processing all coded entries at line " << is.lineNumber()
            << " in file " <<  parentDict.name() << endl;

        ISstream& iss = dynamic_cast<ISstream&>(is);
        std::istream& stdis = iss.stdStream();

        const std::istream::streampos startPos = stdis.tellg();

        // Starting state
        dictionary save(parentDict);

        // Replace entries on the run-time selection table for the
        // function entries
        executeprimitiveEntryIstreamMemberFunctionTable& et =
            *executeprimitiveEntryIstreamMemberFunctionTablePtr_;

        // Preserve old entries
        executeprimitiveEntryIstreamMemberFunctionTable entryTable(et);

        // Replace entries. Keep old ones or not? E.g #include?
        //et.clear();
        forAllConstIter
        (
            preProcesspreProcessPrimitiveEntryMemberFunctionTable,
            *preProcesspreProcessPrimitiveEntryMemberFunctionTablePtr_,
            iter
        )
        {
            //Info<< "Switching to preprocess mode for"
            //    << " primitiveEntry function " << iter.key() << endl;
            et.set(iter.key(), iter());
        }


        executedictionaryIstreamMemberFunctionTable& dt =
            *executedictionaryIstreamMemberFunctionTablePtr_;
        // Preserve old entries
        executedictionaryIstreamMemberFunctionTable dictTable(dt);

        // Replace entries. Keep old ones or not? E.g #include?
        //dt.clear();
        forAllConstIter
        (
            preProcesspreProcessDictionaryMemberFunctionTable,
            *preProcesspreProcessDictionaryMemberFunctionTablePtr_,
            iter
        )
        {
            //Info<< "Switching to preprocess mode for"
            //    << " dictionary function " << iter.key() << endl;
            dt.set(iter.key(), iter());
        }



        // Read dictionary
        parentDict.read(is);


        // Make all dynamicCode. This will now compile in parallel which
        // is the whole reason of this code
        if (Foam::isDir("dynamicCode"))
        {
            system("wmake all dynamicCode");
        }


        // Restore stream

        // Dictionary reading sets bad/fail bit (reads past eof?) so reset
        // state before trying to seek
        stdis.clear();
        // Take over std state to ISstream level
        is.setState(stdis.rdstate());

        stdis.seekg(startPos, stdis.beg);

        // Restore table
        et.transfer(entryTable);
        dt.transfer(dictTable);

        // Restore dictionary
        parentDict.transfer(save);

        //Info<< "#" << typeName
        //    << " pre-processed all coded entres"
        //    << " in file " <<  parentDict.name() << endl;
    }
    return true;
}


// ************************************************************************* //
