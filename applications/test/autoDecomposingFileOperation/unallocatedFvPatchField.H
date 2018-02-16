/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
class unallocatedFvMesh;


// Forward declaration of friend functions and operators

template<class Type>
class unallocatedFvPatchField;

template<class Type>
class calculatedUnallocatedFvPatchField;

template<class Type>
class fvMatrix;

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
        const DimensionedField<Type, unallocatedFvMesh>& internalField_;

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
    typedef calculatedUnallocatedFvPatchField<Type> Calculated;


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
                const DimensionedField<Type, unallocatedFvMesh>& iF
            ),
            (p, iF)
        );

//         declareRunTimeSelectionTable
//         (
//             tmp,
//             unallocatedFvPatchField,
//             patchMapper,
//             (
//                 const unallocatedFvPatchField<Type>& ptf,
//                 const fvPatch& p,
//                 const DimensionedField<Type, unallocatedFvMesh>& iF,
//                 const fvPatchFieldMapper& m
//             ),
//             (dynamic_cast<const unallocatedFvPatchFieldType&>(ptf), p, iF, m)
//         );
        declareRunTimeSelectionTable
        (
            tmp,
            unallocatedFvPatchField,
            patchMapper,
            (
                const fvPatchField<Type>& ptf,
                const fvPatch& p,
                const DimensionedField<Type, unallocatedFvMesh>& iF,
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
                const DimensionedField<Type, unallocatedFvMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&
        );

        //- Construct from patch, internal field and value
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const Type& value
        );

        //- Construct from patch and internal field and patch type
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const word& patchType
        );

        //- Construct from patch and internal field and patch field
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given fvPatchField onto a new
        //  patch
        unallocatedFvPatchField
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
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
            const DimensionedField<Type, unallocatedFvMesh>&
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
            const DimensionedField<Type, unallocatedFvMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given unallocatedFvPatchField mapped onto a new patch
        static tmp<unallocatedFvPatchField<Type>> New
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<unallocatedFvPatchField<Type>> New
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        static tmp<unallocatedFvPatchField<Type>> NewCalculatedType
        (
            const fvPatch&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<unallocatedFvPatchField<Type>> NewCalculatedType
        (
            const unallocatedFvPatchField<Type2>&
        );


    //- Destructor
    virtual ~unallocatedFvPatchField<Type>()
    {}


    // Member functions

//        // Attributes
//
//            //- Return the type of the calculated for of
//            //  unallocatedFvPatchField
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
//            //- Return true if the value of the patch field
//            //  is altered by assignment (the default)
//            virtual bool assignable() const
//            {
//                return true;
//            }
//
//            //- Return true if this patch field is coupled
//            virtual bool coupled() const
//            {
//                return false;
//            }
//
//
//        // Access
//
//            //- Return local objectRegistry
//            const objectRegistry& db() const;

            //- Return patch
            const fvPatch& patch() const
            {
                return patch_;
            }
//
//            //- Return dimensioned internal field reference
//            const DimensionedField<Type, unallocatedFvMesh>&
//            internalField() const
//            {
//                return internalField_;
//            }
//
//            //- Return internal field reference
//            const Field<Type>& primitiveField() const
//            {
//                return internalField_;
//            }
//
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
//
//            //- Return true if the boundary condition has already been updated
//            bool updated() const
//            {
//                return updated_;
//            }
//
//            //- Return true if the matrix has already been manipulated
//            bool manipulatedMatrix() const
//            {
//                return manipulatedMatrix_;
//            }
//
//
//        // Mapping functions
//
//            //- Map (and resize as needed) from self given a mapping object
//            virtual void autoMap
//            (
//                const fvPatchFieldMapper&
//            );
//
//            //- Reverse map the given fvPatchField onto this fvPatchField
//            virtual void rmap
//            (
//                const fvPatchField<Type>&,
//                const labelList&
//            );
//
//
//        // Evaluation functions
//
//            //- Return patch-normal gradient
//            virtual tmp<Field<Type>> snGrad() const;
//
//            //- Return patch-normal gradient for coupled-patches
//            //  using the deltaCoeffs provided
//            virtual tmp<Field<Type>> snGrad
//            (
//                const scalarField& deltaCoeffs
//            ) const
//            {
//                NotImplemented;
//                return *this;
//            }
//
//            //- Update the coefficients associated with the patch field
//            //  Sets Updated to true
//            virtual void updateCoeffs();
//
//            //- Update the coefficients associated with the patch field
//            //  with a weight field (0..1). This weight field is ususally
//            //  provided as the amount of geometric overlap for 'duplicate'
//            //  patches. Sets Updated to true
//            virtual void updateWeightedCoeffs(const scalarField& weights);
//
//            //- Return internal field next to patch as patch field
//            virtual tmp<Field<Type>> patchInternalField() const;
//
//            //- Return internal field next to patch as patch field
//            virtual void patchInternalField(Field<Type>&) const;
//
//            //- Return patchField on the opposite patch of a coupled patch
//            virtual tmp<Field<Type>> patchNeighbourField() const
//            {
//                NotImplemented;
//                return *this;
//            }
//
//            //- Initialise the evaluation of the patch field
//            virtual void initEvaluate
//            (
//                const Pstream::commsTypes commsType =
//                    Pstream::commsTypes::blocking
//            )
//            {}
//
//            //- Evaluate the patch field, sets Updated to false
//            virtual void evaluate
//            (
//                const Pstream::commsTypes commsType =
//                    Pstream::commsTypes::blocking
//            );
//
//
//            //- Return the matrix diagonal coefficients corresponding to the
//            //  evaluation of the value of this patchField with given weights
//            virtual tmp<Field<Type>> valueInternalCoeffs
//            (
//                const tmp<Field<scalar>>&
//            ) const
//            {
//                NotImplemented;
//                return *this;
//            }
//
//            //- Return the matrix source coefficients corresponding to the
//            //  evaluation of the value of this patchField with given weights
//            virtual tmp<Field<Type>> valueBoundaryCoeffs
//            (
//                const tmp<Field<scalar>>&
//            ) const
//            {
//                NotImplemented;
//                return *this;
//            }
//
//            //- Return the matrix diagonal coefficients corresponding to the
//            //  evaluation of the gradient of this patchField
//            virtual tmp<Field<Type>> gradientInternalCoeffs() const
//            {
//                NotImplemented;
//                return *this;
//            }
//
//            //- Return the matrix diagonal coefficients corresponding to the
//            //  evaluation of the gradient of this coupled patchField
//            //  using the deltaCoeffs provided
//            virtual tmp<Field<Type>> gradientInternalCoeffs
//            (
//                const scalarField& deltaCoeffs
//            ) const
//            {
//                NotImplemented;
//                return *this;
//            }
//
//            //- Return the matrix source coefficients corresponding to the
//            //  evaluation of the gradient of this patchField
//            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const
//            {
//                NotImplemented;
//                return *this;
//            }
//
//            //- Return the matrix source coefficients corresponding to the
//            //  evaluation of the gradient of this coupled patchField
//            //  using the deltaCoeffs provided
//            virtual tmp<Field<Type>> gradientBoundaryCoeffs
//            (
//                const scalarField& deltaCoeffs
//            ) const
//            {
//                NotImplemented;
//                return *this;
//            }
//
//
//            //- Manipulate matrix
//            virtual void manipulateMatrix(fvMatrix<Type>& matrix);
//
//            //- Manipulate matrix with given weights
//            virtual void manipulateMatrix
//            (
//                fvMatrix<Type>& matrix,
//                const scalarField& weights
//            );
//
//
        // I-O

            //- Write
            virtual void write(Ostream&) const;
//
//            //- Helper function to write the keyword and entry only if the
//            //  values are not equal. The value is then output as value2
//            template<class EntryType>
//            void writeEntryIfDifferent
//            (
//                Ostream& os,
//                const word& entryName,
//                const EntryType& value1,
//                const EntryType& value2
//            ) const;
//
//
//        // Check
//
//            //- Check fvPatchField<Type> against given fvPatchField<Type>
//            void check(const fvPatchField<Type>&) const;
//
//
//    // Member operators
//
        virtual void operator=(const UList<Type>&);
//
//        virtual void operator=(const fvPatchField<Type>&);
//        virtual void operator+=(const fvPatchField<Type>&);
//        virtual void operator-=(const fvPatchField<Type>&);
//        virtual void operator*=(const fvPatchField<scalar>&);
//        virtual void operator/=(const fvPatchField<scalar>&);
//
//        virtual void operator+=(const Field<Type>&);
//        virtual void operator-=(const Field<Type>&);
//
//        virtual void operator*=(const Field<scalar>&);
//        virtual void operator/=(const Field<scalar>&);
//
//        virtual void operator=(const Type&);
//        virtual void operator+=(const Type&);
//        virtual void operator-=(const Type&);
//        virtual void operator*=(const scalar);
//        virtual void operator/=(const scalar);
//
//
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
    //#include "calculatedFvPatchField.H"
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
