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

Class
    Foam::dynamicCodeContext

Description
    Encapsulation of dynamic code dictionaries. Extracts all strings from
    dictionary, assuming they are C++ code snippets. Constructs a SHA1 key.

SourceFiles
    dynamicCodeContext.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCodeContext_H
#define dynamicCodeContext_H

#include "dictionary.H"
#include "SHA1Digest.H"
#include "OSHA1stream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicCodeContext Declaration
\*---------------------------------------------------------------------------*/

class dynamicCodeContext
{
    // Private data

        //- Collection of filter variables
        HashTable<string> filterVars_;

        //- Running stream
        mutable OSHA1stream stream_;

        //- Calculated SHA1Digest
        mutable autoPtr<SHA1Digest> sha1_;


public:

    // Constructors

        //- Construct from a dictionary. Adds as filterVars
        //  - codeInclude, codeOptions without line directive
        //  - all others strings with line directive
        dynamicCodeContext(const dictionary&);

        //- Construct null
        dynamicCodeContext();


    // Member functions

        //- Collection of filter variables
        const HashTable<string>& filterVars() const
        {
            return filterVars_;
        }

//         //- Helper: return variable content if set, default otherwise
//         const string& filterVar
//         (
//             const word& var,
//             const word& defaultVal = ""
//         ) const;

        //- Return SHA1 digest calculated from current stream. Call only
        //  after all filterVars have been added
        const SHA1Digest& sha1() const
        {
            if (sha1_.empty())
            {
                sha1_.reset(new SHA1Digest(stream_.digest()));
            }
            return sha1_();
        }

        //- Reset the context
        void clear();

        //- Read from dictionary
        void read(const dictionary&);

        //- Helper: write filter variables with original (unexpanded) value
        //  (from dictionary) to stream
        void write(Ostream&, const dictionary&) const;


        //- Helper: add a filter variable from a dictionary. If addDirective
        //          if adds the dictionary line number (after updating the
        //          SHA1)
        void addFilterVariable
        (
            const bool addDirective,
            const dictionary& dict,
            const word& key
        );

        //- Helper: add entries for all variables that occur in (C++) template
        //  file
        void addFilterVariables
        (
            const fileName&,
            const dictionary& dict
        );

        //- Helper: add \#line directive
        static void addLineDirective
        (
            string&,
            const label lineNum,
            const fileName& name
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
