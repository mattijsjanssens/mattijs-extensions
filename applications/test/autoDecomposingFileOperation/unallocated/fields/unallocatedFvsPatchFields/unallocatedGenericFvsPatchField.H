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
    Foam::unallocatedGenericFvsPatchField

Description
    This boundary condition provides a generic version of the \c calculated
    condition, useful as a fallback for handling unknown patch types.  Not
    generally applicable as a user-specified condition.

    It uses parsing logic to decide the type of any Field entry, e.g.
    'nonuniform List<scalar>'. This has problems with 'nonuniform 0' since
    it does not output the type. Normally these get set to be a scalarField
    but this gives problem in reconstruction since there we use in-place
    mapping. What used to be a zero-sized field (can be read a scalar or vector
    or anything) now becomes a non-zero sized scalar field which cannot be
    read by a vectorField (in binary, when it is trying to steal the compound).

See also
    Foam::genericFvPatchField

SourceFiles
    unallocatedGenericFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedGenericFvsPatchField_H
#define unallocatedGenericFvsPatchField_H

#include "HashPtrTable.H"
#include "unallocatedFvsPatchField.H"
#include "fvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class unallocatedGenericFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unallocatedGenericFvsPatchField
:
    public unallocatedFvsPatchField<Type>
{
    // Private data

        word actualTypeName_;
        dictionary dict_;
        bool hasValue_;

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

    //- Runtime type information
    TypeName("generic");


    // Constructors

        //- Construct from patch and internal field
        unallocatedGenericFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedGenericFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        unallocatedGenericFvsPatchField
        (
            const unallocatedFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        unallocatedGenericFvsPatchField
        (
            const unallocatedGenericFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<unallocatedFvsPatchField<Type>> clone() const
        {
            return tmp<unallocatedFvsPatchField<Type>>
            (
                new unallocatedGenericFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unallocatedGenericFvsPatchField
        (
            const unallocatedGenericFvsPatchField<Type>&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<unallocatedFvsPatchField<Type>> clone
        (
            const DimensionedField<Type, unallocatedSurfaceMesh>& iF
        ) const
        {
            return tmp<unallocatedFvsPatchField<Type>>
            (
                new unallocatedGenericFvsPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        const word& actualTypeName() const
        {
            return actualTypeName_;
        }

        //- Extract list of entry names and their types (scalar, vector etc)
        List<Pair<word>> entryTypes() const;

        //- Add dummy field of given type
        void addEntry(const word& key, const word& type);

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        );

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const unallocatedFvsPatchField<Type>&,
            const labelList&
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unallocatedGenericFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
