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
    Foam::unallocatedFvPatchField

Description
    Abstract base class with a fat-interface to all derived classes
    covering all possible ways in which they might be used.

    The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.

    The next level of derivation covers all the specialised types with
    specific evaluation proceedures, particularly with respect to specific
    fields.

SourceFiles
    unallocatedFvPatchField.C
    unallocatedFvPatchFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvPatchField_H
#define unallocatedFvPatchField_H

#include "fvPatch.H"
#include "DimensionedField.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class objectRegistry;
class dictionary;
class fvPatchFieldMapper;
class unallocatedVolMesh;


// Forward declaration of friend functions and operators

template<class Type>
class unallocatedFvPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const unallocatedFvPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                        Class unallocatedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unallocatedFvPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, unallocatedVolMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;

        //- Update index used so that manipulateMatrix is called only once
        //  during the construction of the matrix
        bool manipulatedMatrix_;

        //- Optional patch type, used to allow specified boundary conditions
        //  to be applied to constraint patches by providing the constraint
        //  patch type as 'patchType'
        word patchType_;


public:

    typedef fvPatch Patch;


    //- Runtime type information
    TypeName("unallocatedFvPatchField");

    //- Debug switch to disallow the use of genericunallocatedFvPatchField
    static int disallowGenericunallocatedFvPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            unallocatedFvPatchField,
            patch,
            (
                const fvPatch& p,
                const DimensionedField<Type, unallocatedVolMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            unallocatedFvPatchField,
            patchMapper,
            (
                const unallocatedFvPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, unallocatedVolMesh>& iF,
                const fvPatchFieldMapper& m
            ),
            (ptf, p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            unallocatedFvPatchField,
            dictionary,
            (
                const fvPatch& p,
                const DimensionedField<Type, unallocatedVolMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&
        );

        //- Construct from patch, internal field and value
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const Type& value
        );

        //- Construct from patch and internal field and patch type
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const word& patchType
        );

        //- Construct from patch and internal field and patch field
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given fvPatchField onto a new
        //  patch
        unallocatedFvPatchField
        (
            const unallocatedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        unallocatedFvPatchField(const unallocatedFvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<unallocatedFvPatchField<Type>> clone() const
        {
            return tmp<unallocatedFvPatchField<Type>>
            (
                new unallocatedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unallocatedFvPatchField
        (
            const unallocatedFvPatchField<Type>&,
            const DimensionedField<Type, unallocatedVolMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<unallocatedFvPatchField<Type>> clone
        (
            const DimensionedField<Type, unallocatedVolMesh>& iF
        ) const
        {
            return tmp<unallocatedFvPatchField<Type>>
            (
                new unallocatedFvPatchField<Type>(*this, iF)
            );
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<unallocatedFvPatchField<Type>> New
        (
            const word&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static tmp<unallocatedFvPatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given unallocatedFvPatchField mapped onto a new patch
        static tmp<unallocatedFvPatchField<Type>> New
        (
            const unallocatedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<unallocatedFvPatchField<Type>> New
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const dictionary&
        );


    //- Destructor
    virtual ~unallocatedFvPatchField<Type>()
    {}


    // Member functions

       // Access

            //- Return patch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Optional patch type
            const word& patchType() const
            {
                return patchType_;
            }

            //- Optional patch type
            word& patchType()
            {
                return patchType_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, unallocatedVolMesh>&
            internalField() const
            {
                return internalField_;
            }


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


        // I-O

            //- Write
            virtual void write(Ostream&) const;


   // Member operators

        virtual void operator=(const UList<Type>&);

        // Force an assignment irrespective of form of patch

        virtual void operator==(const unallocatedFvPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const unallocatedFvPatchField<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unallocatedFvPatchField.C"
#endif

/*
#define makeFvPatchField(fvPatchTypeField)                                     \
                                                                               \
defineNamedTemplateTypeNameAndDebug(fvPatchTypeField, 0);                      \
template<>                                                                     \
int fvPatchTypeField::disallowGenericFvPatchField                              \
(                                                                              \
    debug::debugSwitch("disallowGenericFvPatchField", 0)                       \
);                                                                             \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patch);                  \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patchMapper);            \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, dictionary);
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
