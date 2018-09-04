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
    Foam::unallocatedFvsPatchField

Description
    Abstract base class with a fat-interface to all derived classes
    covering all possible ways in which they might be used.

    The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.

    The next level of derivation covers all the specialised types with
    specific evaluation proceedures, particularly with respect to specific
    fields.

SourceFiles
    unallocatedFvsPatchField.C
    unallocatedFvsPatchFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvsPatchField_H
#define unallocatedFvsPatchField_H

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
class unallocatedSurfaceMesh;


// Forward declaration of friend functions and operators

template<class Type>
class unallocatedFvsPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const unallocatedFvsPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                  Class unallocatedFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unallocatedFvsPatchField
:
    public Field<Type>
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, unallocatedSurfaceMesh>& internalField_;


public:

    typedef fvPatch Patch;


    //- Runtime type information
    TypeName("unallocatedFvsPatchField");

//    //- Debug switch to disallow the use of genericFvsPatchField
//    static int disallowGenericUnallocatedFvsPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            unallocatedFvsPatchField,
            patch,
            (
                const fvPatch& p,
                const DimensionedField<Type, unallocatedSurfaceMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            unallocatedFvsPatchField,
            patchMapper,
            (
                const unallocatedFvsPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
                const fvPatchFieldMapper& m
            ),
            (ptf, p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            unallocatedFvsPatchField,
            dictionary,
            (
                const fvPatch& p,
                const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        unallocatedFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&
        );

        //- Construct from patch, internal field and value
        unallocatedFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const Type& value
        );

        //- Construct from patch and internal field and patch type
        unallocatedFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const word& patchType
        );

        //- Construct from patch and internal field and patch field
        unallocatedFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given unallocatedFvsPatchField onto a new
        //  patch
        unallocatedFvsPatchField
        (
            const unallocatedFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        unallocatedFvsPatchField(const unallocatedFvsPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<unallocatedFvsPatchField<Type>> clone() const
        {
            return tmp<unallocatedFvsPatchField<Type>>
            (
                new unallocatedFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unallocatedFvsPatchField
        (
            const unallocatedFvsPatchField<Type>&,
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
                new unallocatedFvsPatchField<Type>(*this, iF)
            );
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<unallocatedFvsPatchField<Type>> New
        (
            const word&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static tmp<unallocatedFvsPatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given unallocatedFvsPatchField mapped onto a new patch
        static tmp<unallocatedFvsPatchField<Type>> New
        (
            const unallocatedFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<unallocatedFvsPatchField<Type>> New
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const dictionary&
        );


    //- Destructor
    virtual ~unallocatedFvsPatchField<Type>()
    {}


    // Member functions

       // Access
//
//            //- Return local objectRegistry
//            const objectRegistry& db() const;

            //- Return patch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, unallocatedSurfaceMesh>&
            internalField() const
            {
                return internalField_;
            }


//            //- Return internal field reference
//            const Field<Type>& primitiveField() const
//            {
//                return internalField_;
//            }
//
//            //- Return the type of the calculated for of fvsPatchField
//            static const word& calculatedType();
//
//            //- Return true if this patch field fixes a value.
//            //  Needed to check if a level has to be specified while solving
//            //  Poissons equations.
//            virtual bool fixesValue() const
//            {
//                return false;
//            }
//
//            //- Return true if this patch field is coupled
//            virtual bool coupled() const
//            {
//                return false;
//            }


        // Mapping functions

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


        // Check

            //- Check unallocatedFvsPatchField<Type> against given
            //  unallocatedFvsPatchField<Type>
            void check(const unallocatedFvsPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const unallocatedFvsPatchField<Type>&);
        virtual void operator+=(const unallocatedFvsPatchField<Type>&);
        virtual void operator-=(const unallocatedFvsPatchField<Type>&);
        virtual void operator*=(const unallocatedFvsPatchField<scalar>&);
        virtual void operator/=(const unallocatedFvsPatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const unallocatedFvsPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const unallocatedFvsPatchField<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unallocatedFvsPatchField.C"
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
