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
    Foam::unallocatedEmptyFvPatchField

Description
    This boundary condition provides a generic version of the \c calculated
    condition, useful as a fallback for handling unknown patch types.  Not
    generally applicable as a user-specified condition.

See also
    Foam::calculatedFvPatchField

SourceFiles
    unallocatedEmptyFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedEmptyFvPatchField_H
#define unallocatedEmptyFvPatchField_H

#include "unallocatedFvPatchField.H"
#include "emptyFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class unallocatedEmptyFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unallocatedEmptyFvPatchField
:
    public unallocatedFvPatchField<Type>
{


public:

    //- Runtime type information
    TypeName(emptyFvPatchField<Type>::typeName_());


    // Constructors

        //- Construct from patch and internal field
        unallocatedEmptyFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedEmptyFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        unallocatedEmptyFvPatchField
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedFvMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        unallocatedEmptyFvPatchField
        (
            const unallocatedEmptyFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<unallocatedFvPatchField<Type>> clone() const
        {
            return tmp<unallocatedFvPatchField<Type>>
            (
                new unallocatedEmptyFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unallocatedEmptyFvPatchField
        (
            const unallocatedEmptyFvPatchField<Type>&,
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
                new unallocatedEmptyFvPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unallocatedEmptyFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
