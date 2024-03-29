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

#include "unallocatedGenericFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "unallocatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Scalar
defineNamedTemplateTypeNameAndDebug(genericUnallocatedFvPatchScalarField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    genericUnallocatedFvPatchScalarField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    genericUnallocatedFvPatchScalarField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    genericUnallocatedFvPatchScalarField,
    patch
);

// Vector
defineNamedTemplateTypeNameAndDebug(genericUnallocatedFvPatchVectorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    genericUnallocatedFvPatchVectorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    genericUnallocatedFvPatchVectorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    genericUnallocatedFvPatchVectorField,
    patch
);

// SphericalTensor
defineNamedTemplateTypeNameAndDebug
(
    genericUnallocatedFvPatchSphericalTensorField,
    0
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    genericUnallocatedFvPatchSphericalTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    genericUnallocatedFvPatchSphericalTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    genericUnallocatedFvPatchSphericalTensorField,
    patch
);

//SymmTensor
defineNamedTemplateTypeNameAndDebug
(
    genericUnallocatedFvPatchSymmTensorField,
    0
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    genericUnallocatedFvPatchSymmTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    genericUnallocatedFvPatchSymmTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    genericUnallocatedFvPatchSymmTensorField,
    patch
);

// Tensor
defineNamedTemplateTypeNameAndDebug(genericUnallocatedFvPatchTensorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    genericUnallocatedFvPatchTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    genericUnallocatedFvPatchTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    genericUnallocatedFvPatchTensorField,
    patch
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
