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

#include "unallocatedProcessorFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "unallocatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Scalar
defineNamedTemplateTypeNameAndDebug(processorUnallocatedFvPatchScalarField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    processorUnallocatedFvPatchScalarField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    processorUnallocatedFvPatchScalarField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchScalarField,
    processorUnallocatedFvPatchScalarField,
    patch
);

// Vector
defineNamedTemplateTypeNameAndDebug(processorUnallocatedFvPatchVectorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    processorUnallocatedFvPatchVectorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    processorUnallocatedFvPatchVectorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchVectorField,
    processorUnallocatedFvPatchVectorField,
    patch
);

// SphericalTensor
defineNamedTemplateTypeNameAndDebug
(
    processorUnallocatedFvPatchSphericalTensorField,
    0
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    processorUnallocatedFvPatchSphericalTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    processorUnallocatedFvPatchSphericalTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSphericalTensorField,
    processorUnallocatedFvPatchSphericalTensorField,
    patch
);

//SymmTensor
defineNamedTemplateTypeNameAndDebug
(
    processorUnallocatedFvPatchSymmTensorField,
    0
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    processorUnallocatedFvPatchSymmTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    processorUnallocatedFvPatchSymmTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchSymmTensorField,
    processorUnallocatedFvPatchSymmTensorField,
    patch
);

// Tensor
defineNamedTemplateTypeNameAndDebug(processorUnallocatedFvPatchTensorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    processorUnallocatedFvPatchTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    processorUnallocatedFvPatchTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvPatchTensorField,
    processorUnallocatedFvPatchTensorField,
    patch
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
