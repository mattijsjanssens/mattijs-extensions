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

#include "unallocatedGenericFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::unallocatedGenericFvPatchField<Type>::unallocatedGenericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF
)
:
    unallocatedFvPatchField<Type>(p, iF),
    genericPatchFieldBase(p.name(), p.size()),
    hasValue_(false)
{
    FatalErrorInFunction
        << "Trying to construct an unallocatedGenericFvPatchField on patch "
        << this->patch().name()
        << abort(FatalError);
}


template<class Type>
Foam::unallocatedGenericFvPatchField<Type>::unallocatedGenericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF,
    const dictionary& dict
)
:
    unallocatedFvPatchField<Type>(p, iF, dict, true),
    genericPatchFieldBase(p.name(), p.size(), dict),
    hasValue_(dict.found("value"))
{}


template<class Type>
Foam::unallocatedGenericFvPatchField<Type>::unallocatedGenericFvPatchField
(
    const unallocatedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    unallocatedFvPatchField<Type>(ptf, p, iF, mapper),
    genericPatchFieldBase
    (
        p.name(),
        p.size(),
        refCast<const unallocatedGenericFvPatchField<Type>>(ptf),
        mapper
    ),
    hasValue_
    (
        refCast<const unallocatedGenericFvPatchField<Type>>(ptf).hasValue_
    )
{}


template<class Type>
Foam::unallocatedGenericFvPatchField<Type>::unallocatedGenericFvPatchField
(
    const unallocatedGenericFvPatchField<Type>& ptf
)
:
    unallocatedFvPatchField<Type>(ptf),
    genericPatchFieldBase(ptf),
    hasValue_(ptf.hasValue_)
{}


template<class Type>
Foam::unallocatedGenericFvPatchField<Type>::unallocatedGenericFvPatchField
(
    const unallocatedGenericFvPatchField<Type>& ptf,
    const DimensionedField<Type, unallocatedVolMesh>& iF
)
:
    unallocatedFvPatchField<Type>(ptf, iF),
    genericPatchFieldBase(ptf),
    hasValue_(ptf.hasValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::unallocatedGenericFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    unallocatedFvPatchField<Type>::autoMap(m);
    genericPatchFieldBase::autoMap(m);
}


template<class Type>
void Foam::unallocatedGenericFvPatchField<Type>::rmap
(
    const unallocatedFvPatchField<Type>& ptf,
    const labelList& addr
)
{
    unallocatedFvPatchField<Type>::rmap(ptf, addr);

    const unallocatedGenericFvPatchField<Type>& dptf =
        refCast<const unallocatedGenericFvPatchField<Type>>(ptf);

    genericPatchFieldBase::rmap(dptf, addr);
}


template<class Type>
void Foam::unallocatedGenericFvPatchField<Type>::write(Ostream& os) const
{
    genericPatchFieldBase::write(os);

    if (hasValue_)
    {
        this->writeEntry("value", os);
    }
}


// ************************************************************************* //
