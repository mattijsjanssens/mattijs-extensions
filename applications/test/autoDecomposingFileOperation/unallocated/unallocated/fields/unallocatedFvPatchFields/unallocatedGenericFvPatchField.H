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

    It uses parsing logic to decide the type of any Field entry, e.g.
    'nonuniform List<scalar>'. This has problems with 'nonuniform 0' since
    it does not output the type. Normally these get set to be a scalarField
    but this gives problem in reconstruction since there we use in-place
    mapping. What used to be a zero-sized field (can be read a scalar or vector
    or anything) now becomes a non-zero sized scalar field which cannot be
    read by a vectorField (in binary, when it is trying to steal the compound).

See also
    Foam::genericPatchFieldBase

SourceFiles
    unallocatedGenericFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedGenericFvPatchField_H
#define unallocatedGenericFvPatchField_H

#include "unallocatedFvPatchField.H"
#include "genericPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class unallocatedGenericFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unallocatedGenericFvPatchField
:
    public unallocatedFvPatchField<Type>,
    public genericPatchFieldBase
{
     // Private data

         const bool hasValue_;


public:

    //- Runtime type information
    TypeName("generic");


    // Constructors

        //- Construct from patch and internal field
        unallocatedGenericFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedGenericFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        unallocatedGenericFvPatchField
        (
            const unallocatedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
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
                new unallocatedGenericFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper&);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const unallocatedFvPatchField<Type>&,
            const labelList&
        );

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
