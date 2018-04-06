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

\*---------------------------------------------------------------------------*/

#include "unallocatedEmptyFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchField.H"
#include "emptyFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::unallocatedEmptyFvPatchField<Type>::unallocatedEmptyFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF
)
:
    unallocatedFvPatchField<Type>(p, iF)
{
    FatalErrorInFunction
        << "Trying to construct an unallocatedEmptyFvPatchField on patch "
        << this->patch().name()
        //<< " of field " << this->internalField().name()
        << abort(FatalError);
}


template<class Type>
Foam::unallocatedEmptyFvPatchField<Type>::unallocatedEmptyFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF,
    const dictionary& dict
)
:
    unallocatedFvPatchField<Type>(p, iF, dict, false)
{}


template<class Type>
Foam::unallocatedEmptyFvPatchField<Type>::unallocatedEmptyFvPatchField
(
    const unallocatedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    unallocatedFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::unallocatedEmptyFvPatchField<Type>::unallocatedEmptyFvPatchField
(
    const unallocatedEmptyFvPatchField<Type>& ptf
)
:
    unallocatedFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::unallocatedEmptyFvPatchField<Type>::unallocatedEmptyFvPatchField
(
    const unallocatedEmptyFvPatchField<Type>& ptf,
    const DimensionedField<Type, unallocatedFvMesh>& iF
)
:
    unallocatedFvPatchField<Type>(ptf, iF)
{}


// ************************************************************************* //
