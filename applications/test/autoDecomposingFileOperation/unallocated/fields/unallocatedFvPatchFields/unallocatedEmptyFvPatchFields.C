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

#include "unallocatedEmptyFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "unallocatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(emptyUnallocatedFvPatchScalarField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    emptyUnallocatedFvPatchScalarField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    emptyUnallocatedFvPatchScalarField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    emptyUnallocatedFvPatchScalarField,
    patch
);

defineNamedTemplateTypeNameAndDebug(emptyUnallocatedFvPatchVectorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    emptyUnallocatedFvPatchVectorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    emptyUnallocatedFvPatchVectorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    emptyUnallocatedFvPatchVectorField,
    patch
);

defineNamedTemplateTypeNameAndDebug
(   emptyUnallocatedFvPatchSphericalTensorField,
    0
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    emptyUnallocatedFvPatchSphericalTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    emptyUnallocatedFvPatchSphericalTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    emptyUnallocatedFvPatchSphericalTensorField,
    patch
);

defineNamedTemplateTypeNameAndDebug(emptyUnallocatedFvPatchSymmTensorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    emptyUnallocatedFvPatchSymmTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    emptyUnallocatedFvPatchSymmTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    emptyUnallocatedFvPatchSymmTensorField,
    patch
);

defineNamedTemplateTypeNameAndDebug(emptyUnallocatedFvPatchTensorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    emptyUnallocatedFvPatchTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    emptyUnallocatedFvPatchTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    emptyUnallocatedFvPatchTensorField,
    patch
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
