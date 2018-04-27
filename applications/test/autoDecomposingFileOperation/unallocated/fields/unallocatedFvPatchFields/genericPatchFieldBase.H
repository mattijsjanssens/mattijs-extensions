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

Class
    Foam::genericPatchFieldBase

Description
    Functionality for 'generic' boundary conditions.

    It uses parsing logic to decide the type of any Field entry, e.g.
    'nonuniform List<scalar>'. This has problems with 'nonuniform 0' since
    it does not output the type. Normally these get set to be a scalarField
    but this gives problem in reconstruction since there we use in-place
    mapping. What used to be a zero-sized field (can be read a scalar or vector
    or anything) now becomes a non-zero sized scalar field which cannot be
    read by a vectorField (in binary, when it is trying to steal the compound).

See also
    Foam::unallocatedGenericFvPatchField
    Foam::unallocatedGenericFvsPatchField

SourceFiles
    genericPatchFieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef genericPatchFieldBase_H
#define genericPatchFieldBase_H

#include "HashPtrTable.H"
#include "primitiveFields.H"
#include "wordPairList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvPatchFieldMapper;

/*---------------------------------------------------------------------------*\
                    Class genericPatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

class genericPatchFieldBase
{
    // Private data

        const word patchName_;
        const label patchSize_;
        word actualTypeName_;
        dictionary dict_;

        HashPtrTable<scalarField> scalarFields_;
        HashPtrTable<vectorField> vectorFields_;
        HashPtrTable<sphericalTensorField> sphericalTensorFields_;
        HashPtrTable<symmTensorField> symmTensorFields_;
        HashPtrTable<tensorField> tensorFields_;


   // Private Member Functions

        void read
        (
            const dictionary& dict,
            HashPtrTable<scalarField>& scalarFields,
            HashPtrTable<vectorField>& vectorFields,
            HashPtrTable<sphericalTensorField>& sphericalTensorFields,
            HashPtrTable<symmTensorField>& symmTensorFields,
            HashPtrTable<tensorField>& tensorFields
        );


public:

    // Constructors

        //- Construct from patch information
        genericPatchFieldBase
        (
            const word& patchName,
            const label patchSize
        );

        //- Construct from patch information and dictionary
        genericPatchFieldBase
        (
            const word& patchName,
            const label patchSize,
            const dictionary&
        );

        //- Construct by mapping given genericPatchFieldBase onto a new patch
        genericPatchFieldBase
        (
            const word& patchName,
            const label patchSize,
            const genericPatchFieldBase&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        genericPatchFieldBase(const genericPatchFieldBase&);


    // Member functions

        const word& actualTypeName() const
        {
            return actualTypeName_;
        }


        // Handling types of generic fields

            //- Extract list of entry names and their types (scalar, vector etc)
            wordPairList entryTypes() const;

            //- Add dummy field of given type
            void addEntry(const word& key, const word& type);

            //- Find matching entry
            static label findEntry
            (
                const wordPairList& lst,
                const word& key
            );

            //- Merge into table from typename to entry list
            static void mergeEntries
            (
                const word& actualTypeName,
                const wordPairList& entryTypes,
                HashTable<wordPairList>& patchFieldToTypes
            );


        //- Map (and resize as needed) from self given a mapping object
        void autoMap
        (
            const fvPatchFieldMapper&
        );

        //- Reverse map the given fvPatchField onto this fvPatchField
        void rmap
        (
            const genericPatchFieldBase&,
            const labelList&
        );

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
