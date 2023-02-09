/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::dummyPrimitiveEntry

Description
    primitiveEntry that prints contents. Used to insert e.g. comment
    into output

SourceFiles
    dummyPrimitiveEntry.C
    dummyPrimitiveEntryIO.C

SeeAlso
    commentDirective.H

\*---------------------------------------------------------------------------*/

#ifndef dummyPrimitiveEntry_H
#define dummyPrimitiveEntry_H

#include "primitiveEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
//class dictionary;

/*---------------------------------------------------------------------------*\
                       Class dummyPrimitiveEntry Declaration
\*---------------------------------------------------------------------------*/

class dummyPrimitiveEntry
:
    public primitiveEntry
{
    // Private Member Functions


public:

    // Constructors

//        //- Construct from keyword and no tokens.
//        //  Contents to be filled with a later assignment
//        explicit dummyPrimitiveEntry(const keyType& key);

        //- Construct from keyword and a single token
        dummyPrimitiveEntry(const keyType& key, const token& tok);

//        //- Construct from keyword and a list of tokens
//        dummyPrimitiveEntry(const keyType& key, const UList<token>& tokens);
//
//        //- Construct from keyword and by transferring a list of tokens
//        dummyPrimitiveEntry(const keyType& key, List<token>&& tokens);
//
//        //- Construct from keyword and ITstream tokens
//        dummyPrimitiveEntry(const keyType& key, const ITstream& is);
//
//        //- Construct from keyword and Istream
//        dummyPrimitiveEntry(const keyType& key, Istream& is);
//
//        //- Construct from keyword, parent dictionary and Istream
//        dummyPrimitiveEntry(const keyType& key, const dictionary& dict, Istream& is);

        //- Construct from keyword and a value. Uses string stream serialization
        template<class T>
        dummyPrimitiveEntry(const keyType& key, const T& val);

//        //- Clone the entry
//        autoPtr<entry> clone(const dictionary&) const
//        {
//            return autoPtr<entry>(new dummyPrimitiveEntry(*this));
//        }


    // Member Functions


        //- Write
        virtual void write(Ostream& os) const
        {
            ITstream& tokens = stream();
            for (const auto& t : tokens)
            {
                //os.writeQuoted("//", false);
                os.writeQuoted(t.stringToken(), false);
                os.write('\n');
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dummyPrimitiveEntryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
