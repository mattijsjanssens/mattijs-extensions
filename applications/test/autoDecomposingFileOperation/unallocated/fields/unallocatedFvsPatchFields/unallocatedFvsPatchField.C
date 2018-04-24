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

#include "dictionary.H"
#include "fvPatchFieldMapper.H"
#include "unallocatedFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::unallocatedFvsPatchField<Type>::unallocatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
Foam::unallocatedFvsPatchField<Type>::unallocatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
Foam::unallocatedFvsPatchField<Type>::unallocatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
    const dictionary& dict
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF)
{
    if (dict.found("value"))
    {
        unallocatedFvsPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        //IOWarningInFunction
        //(
        //    dict
        //)   << "Essential entry 'value' missing."
        //    << " Initialising to zero instead"
        //    << endl;
        Field<Type>::operator=(Zero);
    }
}


template<class Type>
Foam::unallocatedFvsPatchField<Type>::unallocatedFvsPatchField
(
    const unallocatedFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    Field<Type>(ptf, mapper),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
Foam::unallocatedFvsPatchField<Type>::unallocatedFvsPatchField
(
    const unallocatedFvsPatchField<Type>& ptf
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_)
{}


template<class Type>
Foam::unallocatedFvsPatchField<Type>::unallocatedFvsPatchField
(
    const unallocatedFvsPatchField<Type>& ptf,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class Type>
//const Foam::objectRegistry& Foam::unallocatedFvsPatchField<Type>::db() const
//{
//    return patch_.boundaryMesh().mesh();
//}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::check
(
    const unallocatedFvsPatchField<Type>& ptf
) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorInFunction
            << "different patches for unallocatedFvsPatchField<Type>s"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::rmap
(
    const unallocatedFvsPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator=
(
    const unallocatedFvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator+=
(
    const unallocatedFvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator-=
(
    const unallocatedFvsPatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator*=
(
    const unallocatedFvsPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorInFunction
            << "incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator*=(ptf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator/=
(
    const unallocatedFvsPatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    Field<Type>::operator/=(ptf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator==
(
    const unallocatedFvsPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void Foam::unallocatedFvsPatchField<Type>::operator==
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
    const unallocatedFvsPatchField<Type>& ptf
)
{
    ptf.write(os);

    os.check
    (
        "Ostream& operator<<(Ostream&, const unallocatedFvsPatchField<Type>&"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "unallocatedFvsPatchFieldNew.C"

// ************************************************************************* //
