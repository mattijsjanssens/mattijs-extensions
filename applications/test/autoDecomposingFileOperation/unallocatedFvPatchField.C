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

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "dictionary.H"
//#include "fvMesh.H"
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
    const unallocatedFvPatchField<Type>& ptf,
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
    patchType_(ptf.patchType_)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        unallocatedFvPatchField<Type>::operator=(this->patchInternalField());
    }
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

//template<class Type>
//const Foam::objectRegistry& Foam::unallocatedFvPatchField<Type>::db() const
//{
//    return patch_.boundaryMesh().mesh();
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::check(const unallocatedFvPatchField<Type>& ptf) const
//{
//    if (&patch_ != &(ptf.patch_))
//    {
//        FatalErrorInFunction
//            << "different patches for unallocatedFvPatchField<Type>s"
//            << abort(FatalError);
//    }
//}
//
//
//template<class Type>
//Foam::tmp<Foam::Field<Type>> Foam::unallocatedFvPatchField<Type>::snGrad() const
//{
//    return patch_.deltaCoeffs()*(*this - patchInternalField());
//}
//
//
//template<class Type>
//Foam::tmp<Foam::Field<Type>>
//Foam::unallocatedFvPatchField<Type>::patchInternalField() const
//{
//    return patch_.patchInternalField(internalField_);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::patchInternalField(Field<Type>& pif) const
//{
//    patch_.patchInternalField(internalField_, pif);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::autoMap
//(
//    const unallocatedFvPatchFieldMapper& mapper
//)
//{
//    Field<Type>& f = *this;
//
//    if (!this->size() && !mapper.distributed())
//    {
//        f.setSize(mapper.size());
//        if (f.size())
//        {
//            f = this->patchInternalField();
//        }
//    }
//    else
//    {
//        // Map all faces provided with mapping data
//        Field<Type>::autoMap(mapper);
//
//        // For unmapped faces set to internal field value (zero-gradient)
//        if (mapper.hasUnmapped())
//        {
//            Field<Type> pif(this->patchInternalField());
//
//            if
//            (
//                mapper.direct()
//             && notNull(mapper.directAddressing())
//             && mapper.directAddressing().size()
//            )
//            {
//                const labelList& mapAddressing = mapper.directAddressing();
//
//                forAll(mapAddressing, i)
//                {
//                    if (mapAddressing[i] < 0)
//                    {
//                        f[i] = pif[i];
//                    }
//                }
//            }
//            else if (!mapper.direct() && mapper.addressing().size())
//            {
//                const labelListList& mapAddressing = mapper.addressing();
//
//                forAll(mapAddressing, i)
//                {
//                    const labelList& localAddrs = mapAddressing[i];
//
//                    if (!localAddrs.size())
//                    {
//                        f[i] = pif[i];
//                    }
//                }
//            }
//        }
//    }
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::rmap
//(
//    const unallocatedFvPatchField<Type>& ptf,
//    const labelList& addr
//)
//{
//    Field<Type>::rmap(ptf, addr);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::updateCoeffs()
//{
//    updated_ = true;
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::updateWeightedCoeffs(const scalarField& weights)
//{
//    // Default behaviour ignores the weights
//    if (!updated_)
//    {
//        updateCoeffs();
//
//        updated_ = true;
//    }
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
//{
//    if (!updated_)
//    {
//        updateCoeffs();
//    }
//
//    updated_ = false;
//    manipulatedMatrix_ = false;
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::manipulateMatrix(fvMatrix<Type>& matrix)
//{
//    manipulatedMatrix_ = true;
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::manipulateMatrix
//(
//    fvMatrix<Type>& matrix,
//    const scalarField& weights
//)
//{
//    manipulatedMatrix_ = true;
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::write(Ostream& os) const
//{
//    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
//
//    if (patchType_.size())
//    {
//        os.writeKeyword("patchType") << patchType_
//            << token::END_STATEMENT << nl;
//    }
//}
//
//
//template<class Type>
//template<class EntryType>
//void Foam::unallocatedFvPatchField<Type>::writeEntryIfDifferent
//(
//    Ostream& os,
//    const word& entryName,
//    const EntryType& value1,
//    const EntryType& value2
//) const
//{
//    if (value1 != value2)
//    {
//        os.writeKeyword(entryName) << value2 << token::END_STATEMENT << nl;
//    }
//}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator=
//(
//    const UList<Type>& ul
//)
//{
//    Field<Type>::operator=(ul);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator=
//(
//    const unallocatedFvPatchField<Type>& ptf
//)
//{
//    check(ptf);
//    Field<Type>::operator=(ptf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator+=
//(
//    const unallocatedFvPatchField<Type>& ptf
//)
//{
//    check(ptf);
//    Field<Type>::operator+=(ptf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator-=
//(
//    const unallocatedFvPatchField<Type>& ptf
//)
//{
//    check(ptf);
//    Field<Type>::operator-=(ptf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator*=
//(
//    const unallocatedFvPatchField<scalar>& ptf
//)
//{
//    if (&patch_ != &ptf.patch())
//    {
//        FatalErrorInFunction
//            << "incompatible patches for patch fields"
//            << abort(FatalError);
//    }
//
//    Field<Type>::operator*=(ptf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator/=
//(
//    const unallocatedFvPatchField<scalar>& ptf
//)
//{
//    if (&patch_ != &ptf.patch())
//    {
//        FatalErrorInFunction
//            << abort(FatalError);
//    }
//
//    Field<Type>::operator/=(ptf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator+=
//(
//    const Field<Type>& tf
//)
//{
//    Field<Type>::operator+=(tf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator-=
//(
//    const Field<Type>& tf
//)
//{
//    Field<Type>::operator-=(tf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator*=
//(
//    const scalarField& tf
//)
//{
//    Field<Type>::operator*=(tf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator/=
//(
//    const scalarField& tf
//)
//{
//    Field<Type>::operator/=(tf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator=
//(
//    const Type& t
//)
//{
//    Field<Type>::operator=(t);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator+=
//(
//    const Type& t
//)
//{
//    Field<Type>::operator+=(t);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator-=
//(
//    const Type& t
//)
//{
//    Field<Type>::operator-=(t);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator*=
//(
//    const scalar s
//)
//{
//    Field<Type>::operator*=(s);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator/=
//(
//    const scalar s
//)
//{
//    Field<Type>::operator/=(s);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator==
//(
//    const unallocatedFvPatchField<Type>& ptf
//)
//{
//    Field<Type>::operator=(ptf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator==
//(
//    const Field<Type>& tf
//)
//{
//    Field<Type>::operator=(tf);
//}
//
//
//template<class Type>
//void Foam::unallocatedFvPatchField<Type>::operator==
//(
//    const Type& t
//)
//{
//    Field<Type>::operator=(t);
//}
//
//
//// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //
//
//template<class Type>
//Foam::Ostream& Foam::operator<<(Ostream& os, const unallocatedFvPatchField<Type>& ptf)
//{
//    ptf.write(os);
//
//    os.check("Ostream& operator<<(Ostream&, const unallocatedFvPatchField<Type>&");
//
//    return os;
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "unallocatedFvPatchFieldNew.C"

// ************************************************************************* //
