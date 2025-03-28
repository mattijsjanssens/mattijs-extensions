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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::unallocatedFvPatchField<Type>>
Foam::unallocatedFvPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF
)
{
    if (!patchConstructorTablePtr_)
    {
        FatalErrorInFunction << "No constructors-from-patch available"
            << " when constructing from patch " << p.name()
            << " from field type " << patchFieldType
            << exit(FatalError);
    }
    if (debug)
    {
        InfoInFunction
            << "patchFieldType = " << patchFieldType
            << " : " << p.type()
            << endl;
    }

    typename patchConstructorTable::iterator cstrIter =
        patchConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == patchConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown patchField type "
            << patchFieldType << nl << nl
            << "Valid patchField types are :" << endl
            << patchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTablePtr_->find(p.type());

    if
    (
        actualPatchType == word::null
     || actualPatchType != p.type()
    )
    {
        if (patchTypeCstrIter != patchConstructorTablePtr_->end())
        {
            return patchTypeCstrIter()(p, iF);
        }
        else
        {
            return cstrIter()(p, iF);
        }
    }
    else
    {
        tmp<unallocatedFvPatchField<Type>> tfvp = cstrIter()(p, iF);

        // Check if constraint type override and store patchType if so
        if ((patchTypeCstrIter != patchConstructorTablePtr_->end()))
        {
            tfvp.ref().patchType() = actualPatchType;
        }
        return tfvp;
    }
}


template<class Type>
Foam::tmp<Foam::unallocatedFvPatchField<Type>>
Foam::unallocatedFvPatchField<Type>::New
(
    const word& patchFieldType,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::tmp<Foam::unallocatedFvPatchField<Type>>
Foam::unallocatedFvPatchField<Type>::New
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF,
    const dictionary& dict
)
{
    const word patchFieldType(dict.lookup("type"));

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorInFunction << "No constructors-from-dictionary available"
            << " when constructing for patch " << p.name()
            << " from dictionary " << dict
            << exit(FatalError);
    }

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        cstrIter = dictionaryConstructorTablePtr_->find("generic");

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << "Valid patchField types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    if
    (
       !dict.found("patchType")
     || word(dict.lookup("patchType")) != p.type()
    )
    {
        typename dictionaryConstructorTable::iterator patchTypeCstrIter
            = dictionaryConstructorTablePtr_->find(p.type());

        if
        (
            patchTypeCstrIter != dictionaryConstructorTablePtr_->end()
         && patchTypeCstrIter() != cstrIter()
        )
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "inconsistent patch and patchField types for \n"
                   "    patch type " << p.type()
                << " and patchField type " << patchFieldType
                << exit(FatalIOError);
        }
    }

    return cstrIter()(p, iF, dict);
}


template<class Type>
Foam::tmp<Foam::unallocatedFvPatchField<Type>>
Foam::unallocatedFvPatchField<Type>::New
(
    const unallocatedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF,
    const fvPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        InfoInFunction << "Constructing unallocatedFvPatchField<Type>" << endl;
    }

    if (!patchMapperConstructorTablePtr_)
    {
        FatalErrorInFunction << "No constructors-from-mapper available"
            << " when constructing for patch " << p.name()
            << " from fvPatchField " << ptf.type()
            << exit(FatalError);
    }

    typename patchMapperConstructorTable::iterator cstrIter =
        patchMapperConstructorTablePtr_->find(ptf.type());

    if (cstrIter == patchMapperConstructorTablePtr_->end())
    {
        cstrIter = patchMapperConstructorTablePtr_->find("generic");

        if (cstrIter == patchMapperConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown patchField type " << ptf.type() << nl << nl
                << "Valid patchField types are :" << endl
                << patchMapperConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
    }

    return cstrIter()(ptf, p, iF, pfMapper);
}


// ************************************************************************* //
