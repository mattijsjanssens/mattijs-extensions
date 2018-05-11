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

#include "unallocatedProcessorFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::unallocatedProcessorFvPatchField<Type>::unallocatedProcessorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF
)
:
    unallocatedFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::unallocatedProcessorFvPatchField<Type>::unallocatedProcessorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF,
    const dictionary& dict
)
:
    unallocatedFvPatchField<Type>(p, iF, dict, true)
{}


template<class Type>
Foam::unallocatedProcessorFvPatchField<Type>::unallocatedProcessorFvPatchField
(
    const unallocatedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, unallocatedVolMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    unallocatedFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::unallocatedProcessorFvPatchField<Type>::unallocatedProcessorFvPatchField
(
    const unallocatedProcessorFvPatchField<Type>& ptf
)
:
    unallocatedFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::unallocatedProcessorFvPatchField<Type>::unallocatedProcessorFvPatchField
(
    const unallocatedProcessorFvPatchField<Type>& ptf,
    const DimensionedField<Type, unallocatedVolMesh>& iF
)
:
    unallocatedFvPatchField<Type>(ptf, iF)
{}


// ************************************************************************* //
