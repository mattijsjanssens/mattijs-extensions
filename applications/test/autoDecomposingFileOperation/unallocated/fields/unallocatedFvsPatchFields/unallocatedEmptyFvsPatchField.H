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
    Foam::unallocatedEmptyFvsPatchField

Description
    This boundary condition provides a generic version of the \c calculated
    condition, useful as a fallback for handling unknown patch types.  Not
    generally applicable as a user-specified condition.

See also
    Foam::calculatedFvsPatchField

SourceFiles
    unallocatedEmptyFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedEmptyFvsPatchField_H
#define unallocatedEmptyFvsPatchField_H

#include "unallocatedFvsPatchField.H"
#include "emptyFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class unallocatedEmptyFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unallocatedEmptyFvsPatchField
:
    public unallocatedFvsPatchField<Type>
{
public:

    //- Runtime type information
    TypeName(emptyFvsPatchField<Type>::typeName_());


    // Constructors

        //- Construct from patch and internal field
        unallocatedEmptyFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedEmptyFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        unallocatedEmptyFvsPatchField
        (
            const unallocatedFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedSurfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        unallocatedEmptyFvsPatchField
        (
            const unallocatedEmptyFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<unallocatedFvsPatchField<Type>> clone() const
        {
            return tmp<unallocatedFvsPatchField<Type>>
            (
                new unallocatedEmptyFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unallocatedEmptyFvsPatchField
        (
            const unallocatedEmptyFvsPatchField<Type>&,
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
                new unallocatedEmptyFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unallocatedEmptyFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
