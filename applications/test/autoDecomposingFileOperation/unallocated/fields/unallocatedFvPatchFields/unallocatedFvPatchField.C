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

#include "dictionary.H"
#include "fvPatchFieldMapper.H"
#include "unallocatedFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::unallocatedFvPatchField<Type>::unallocatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null)
{}


template<class Type>
Foam::unallocatedFvPatchField<Type>::unallocatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF,
    const Type& value
)
:
    Field<Type>(p.size(), value),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null)
{}


template<class Type>
Foam::unallocatedFvPatchField<Type>::unallocatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF,
    const word& patchType
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(patchType)
{}


template<class Type>
Foam::unallocatedFvPatchField<Type>::unallocatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(word::null)
{}


template<class Type>
Foam::unallocatedFvPatchField<Type>::unallocatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(dict.lookupOrDefault<word>("patchType", word::null))
{
    if (valueRequired)
    {
        if (dict.found("value"))
        {
            Field<Type>::operator=
            (
                Field<Type>("value", dict, p.size())
            );
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Essential entry 'value' missing"
                << exit(FatalIOError);
        }
    }
}


template<class Type>
Foam::unallocatedFvPatchField<Type>::unallocatedFvPatchField
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedFvMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.type())
{
    // Map the value field. Note that we cannot access the patch addressing
    // to do e.g. patchInternalField for initialisation of unmapped values
    this->map(ptf, mapper);
}


template<class Type>
Foam::unallocatedFvPatchField<Type>::unallocatedFvPatchField
(
    const unallocatedFvPatchField<Type>& ptf
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_)
{}


template<class Type>
Foam::unallocatedFvPatchField<Type>::unallocatedFvPatchField
(
    const unallocatedFvPatchField<Type>& ptf,
    const DimensionedField<Type, unallocatedFvMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF),
    updated_(false),
    manipulatedMatrix_(false),
    patchType_(ptf.patchType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::unallocatedFvPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;

    if (patchType_.size() && patchType_ != this->type())
    {
        os.writeKeyword("patchType") << patchType_
            << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::unallocatedFvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::unallocatedFvPatchField<Type>::operator==
(
    const unallocatedFvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::unallocatedFvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::unallocatedFvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const unallocatedFvPatchField<Type>& ptf
)
{
    ptf.write(os);

    os.check
    (
        "Ostream& operator<<(Ostream&, const unallocatedFvPatchField<Type>&"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "unallocatedFvPatchFieldNew.C"

// ************************************************************************* //
