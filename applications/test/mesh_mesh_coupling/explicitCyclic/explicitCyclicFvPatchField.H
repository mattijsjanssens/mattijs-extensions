/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::explicitCyclicFvPatchField

Description
    This boundary condition enforces a cyclic condition between a pair of
    boundaries.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            cyclic;
    }
    \endverbatim

Note
    The patches must be topologically similar, i.e. if the owner patch is
    transformed to the neighbour patch, the patches should be identical (or
    very similar).

SourceFiles
    explicitCyclicFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef explicitCyclicFvPatchField_H
#define explicitCyclicFvPatchField_H

#include "coupledFvPatchField.H"
#include "cyclicLduInterfaceField.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class explicitCyclicFvPatchField
:
    public coupledFvPatchField<Type>,
    public cyclicLduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicFvPatch& cyclicPatch_;


public:

    //- Runtime type information
    TypeName("explicitCyclic");


    // Constructors

        //- Construct from patch and internal field
        explicitCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        explicitCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given explicitCyclicFvPatchField onto a new patch
        explicitCyclicFvPatchField
        (
            const explicitCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        explicitCyclicFvPatchField
        (
            const explicitCyclicFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new explicitCyclicFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        explicitCyclicFvPatchField
        (
            const explicitCyclicFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new explicitCyclicFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const cyclicFvPatch& cyclicPatch() const
            {
                return cyclicPatch_;
            }


        // Evaluation functions

            //- Return neighbour coupled internal cell data
            tmp<Field<Type>> patchNeighbourField() const;

            //- Return reference to neighbour patchField
            const explicitCyclicFvPatchField<Type>& neighbourPatchField() const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& result,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;


        // Cyclic coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return !(cyclicPatch_.parallel() || pTraits<Type>::rank == 0);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicPatch_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "explicitCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
