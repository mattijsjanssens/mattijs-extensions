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

\*---------------------------------------------------------------------------*/

#include "unallocatedGenericFvsPatchField.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::unallocatedGenericFvsPatchField<Type>::unallocatedGenericFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF
)
:
    unallocatedFvsPatchField<Type>(p, iF),
    genericPatchFieldBase(p.name(), p.size()),
    hasValue_(false)
{}


template<class Type>
Foam::unallocatedGenericFvsPatchField<Type>::unallocatedGenericFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
    const dictionary& dict
)
:
    unallocatedFvsPatchField<Type>(p, iF, dict),
    genericPatchFieldBase(p.name(), p.size(), dict),
    hasValue_(dict.found("value"))
{}


template<class Type>
Foam::unallocatedGenericFvsPatchField<Type>::unallocatedGenericFvsPatchField
(
    const unallocatedFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    unallocatedFvsPatchField<Type>(ptf, p, iF, mapper),
    genericPatchFieldBase
    (
        p.name(),
        p.size(),
        refCast<const unallocatedGenericFvsPatchField<Type>>(ptf),
        mapper
    ),
    hasValue_
    (
        refCast<const unallocatedGenericFvsPatchField<Type>>(ptf).hasValue_
    )
{}


template<class Type>
Foam::unallocatedGenericFvsPatchField<Type>::unallocatedGenericFvsPatchField
(
    const unallocatedGenericFvsPatchField<Type>& ptf
)
:
    unallocatedFvsPatchField<Type>(ptf),
    genericPatchFieldBase(ptf),
    hasValue_(ptf.hasValue_)
{}


template<class Type>
Foam::unallocatedGenericFvsPatchField<Type>::unallocatedGenericFvsPatchField
(
    const unallocatedGenericFvsPatchField<Type>& ptf,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF
)
:
    unallocatedFvsPatchField<Type>(ptf, iF),
    genericPatchFieldBase(ptf),
    hasValue_(ptf.hasValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::unallocatedGenericFvsPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    unallocatedFvsPatchField<Type>::autoMap(m);
    genericPatchFieldBase::autoMap(m);
}


template<class Type>
void Foam::unallocatedGenericFvsPatchField<Type>::rmap
(
    const unallocatedFvsPatchField<Type>& ptf,
    const labelList& addr
)
{
    unallocatedFvsPatchField<Type>::rmap(ptf, addr);

    const unallocatedGenericFvsPatchField<Type>& dptf =
        refCast<const unallocatedGenericFvsPatchField<Type>>(ptf);

    genericPatchFieldBase::rmap(dptf, addr);
}


template<class Type>
void Foam::unallocatedGenericFvsPatchField<Type>::write(Ostream& os) const
{
    genericPatchFieldBase::write(os);

    if (hasValue_)
    {
        this->writeEntry("value", os);
    }
}


// ************************************************************************* //
