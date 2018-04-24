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

#include "unallocatedEmptyFvsPatchField.H"
#include "fvPatchFieldMapper.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::unallocatedEmptyFvsPatchField<Type>::unallocatedEmptyFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF
)
:
    unallocatedFvsPatchField<Type>(p, iF, Field<Type>(0))
{}


template<class Type>
Foam::unallocatedEmptyFvsPatchField<Type>::unallocatedEmptyFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
    const dictionary& dict
)
:
    unallocatedFvsPatchField<Type>(p, iF, Field<Type>(0))
{
//    if (!isType<emptyFvPatch>(p))
//    {
//        FatalIOErrorInFunction
//        (
//            dict
//        )   << "patch " << p.name() << " not empty type. "
//            << "Patch type = " << p.type()
//            << exit(FatalIOError);
//    }
}


template<class Type>
Foam::unallocatedEmptyFvsPatchField<Type>::unallocatedEmptyFvsPatchField
(
    const unallocatedFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    unallocatedFvsPatchField<Type>(p, iF, Field<Type>(0))
{
//    if (!isType<emptyFvPatch>(this->patch()))
//    {
//        FatalErrorInFunction
//            << "Field type does not correspond to patch type for patch "
//            << this->patch().index() << "." << endl
//            << "Field type: " << typeName << endl
//            << "Patch type: " << this->patch().type()
//            << exit(FatalError);
//    }
}


template<class Type>
Foam::unallocatedEmptyFvsPatchField<Type>::unallocatedEmptyFvsPatchField
(
    const unallocatedEmptyFvsPatchField<Type>& ptf
)
:
    unallocatedFvsPatchField<Type>
    (
        ptf.patch(),
        ptf.internalField(),
        Field<Type>(0)
    )
{}


template<class Type>
Foam::unallocatedEmptyFvsPatchField<Type>::unallocatedEmptyFvsPatchField
(
    const unallocatedEmptyFvsPatchField<Type>& ptf,
    const DimensionedField<Type, unallocatedSurfaceMesh>& iF
)
:
    unallocatedFvsPatchField<Type>(ptf.patch(), iF, Field<Type>(0))
{}


// ************************************************************************* //
