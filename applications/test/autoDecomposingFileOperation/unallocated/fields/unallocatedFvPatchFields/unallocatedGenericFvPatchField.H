/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::unallocatedGenericFvPatchField

Description
    This boundary condition provides a generic version of the \c calculated
    condition, useful as a fallback for handling unknown patch types.  Not
    generally applicable as a user-specified condition.

See also
    Foam::genericFvPatchField

SourceFiles
    unallocatedGenericFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedGenericFvPatchField_H
#define unallocatedGenericFvPatchField_H

#include "HashPtrTable.H"
#include "unallocatedFvPatchField.H"
#include "unallocatedFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class genericFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unallocatedGenericFvPatchField
:
    public unallocatedFvPatchField<Type>
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
        unallocatedGenericFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedGenericFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        unallocatedGenericFvPatchField
        (
            const unallocatedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        unallocatedGenericFvPatchField
        (
            const unallocatedGenericFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<unallocatedFvPatchField<Type>> clone() const
        {
            return tmp<unallocatedFvPatchField<Type>>
            (
                new unallocatedGenericFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unallocatedGenericFvPatchField
        (
            const unallocatedGenericFvPatchField<Type>&,
            const DimensionedField<Type, unallocatedFvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<unallocatedFvPatchField<Type>> clone
        (
            const DimensionedField<Type, unallocatedFvMesh>& iF
        ) const
        {
            return tmp<unallocatedFvPatchField<Type>>
            (
                new unallocatedGenericFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const unallocatedFvPatchField<Type>&,
                const labelList&
            );



        // Helper: select storage

            static HashPtrTable<Field<scalar>>& fieldTable
            (
                HashPtrTable<scalarField>& scalarFields,
                HashPtrTable<vectorField>& vectorFields,
                HashPtrTable<sphericalTensorField>& sphericalTensorFields,
                HashPtrTable<symmTensorField>& symmTensorFields,
                HashPtrTable<tensorField>& tensorFields,
                const scalar&
            )
            {
                return scalarFields;
            }
            static HashPtrTable<Field<vector>>& fieldTable
            (
                HashPtrTable<scalarField>& scalarFields,
                HashPtrTable<vectorField>& vectorFields,
                HashPtrTable<sphericalTensorField>& sphericalTensorFields,
                HashPtrTable<symmTensorField>& symmTensorFields,
                HashPtrTable<tensorField>& tensorFields,
                const vector&
            )
            {
                return vectorFields;
            }
            static HashPtrTable<Field<sphericalTensor>>& fieldTable
            (
                HashPtrTable<scalarField>& scalarFields,
                HashPtrTable<vectorField>& vectorFields,
                HashPtrTable<sphericalTensorField>& sphericalTensorFields,
                HashPtrTable<symmTensorField>& symmTensorFields,
                HashPtrTable<tensorField>& tensorFields,
                const sphericalTensor&
            )
            {
                return sphericalTensorFields;
            }
            static HashPtrTable<Field<symmTensor>>& fieldTable
            (
                HashPtrTable<scalarField>& scalarFields,
                HashPtrTable<vectorField>& vectorFields,
                HashPtrTable<sphericalTensorField>& sphericalTensorFields,
                HashPtrTable<symmTensorField>& symmTensorFields,
                HashPtrTable<tensorField>& tensorFields,
                const symmTensor&
            )
            {
                return symmTensorFields;
            }
            static HashPtrTable<Field<tensor>>& fieldTable
            (
                HashPtrTable<scalarField>& scalarFields,
                HashPtrTable<vectorField>& vectorFields,
                HashPtrTable<sphericalTensorField>& sphericalTensorFields,
                HashPtrTable<symmTensorField>& symmTensorFields,
                HashPtrTable<tensorField>& tensorFields,
                const tensor&
            )
            {
                return tensorFields;
            }
            template<class FldType>
            HashPtrTable<Field<FldType>>& fieldTable()
            {
                return fieldTable
                (
                    scalarFields_,
                    vectorFields_,
                    sphericalTensorFields_,
                    symmTensorFields_,
                    tensorFields_,
                    pTraits<FldType>::zero
                );
            }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unallocatedGenericFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
