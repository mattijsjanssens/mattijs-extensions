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

#include "unallocatedProcessorFvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "unallocatedFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Scalar
defineNamedTemplateTypeNameAndDebug(processorUnallocatedFvsPatchScalarField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchScalarField,
    processorUnallocatedFvsPatchScalarField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchScalarField,
    processorUnallocatedFvsPatchScalarField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchScalarField,
    processorUnallocatedFvsPatchScalarField,
    patch
);

// Vector
defineNamedTemplateTypeNameAndDebug(processorUnallocatedFvsPatchVectorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchVectorField,
    processorUnallocatedFvsPatchVectorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchVectorField,
    processorUnallocatedFvsPatchVectorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchVectorField,
    processorUnallocatedFvsPatchVectorField,
    patch
);

// SphericalTensor
defineNamedTemplateTypeNameAndDebug
(
    processorUnallocatedFvsPatchSphericalTensorField,
    0
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchSphericalTensorField,
    processorUnallocatedFvsPatchSphericalTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchSphericalTensorField,
    processorUnallocatedFvsPatchSphericalTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchSphericalTensorField,
    processorUnallocatedFvsPatchSphericalTensorField,
    patch
);

//SymmTensor
defineNamedTemplateTypeNameAndDebug
(
    processorUnallocatedFvsPatchSymmTensorField,
    0
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchSymmTensorField,
    processorUnallocatedFvsPatchSymmTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchSymmTensorField,
    processorUnallocatedFvsPatchSymmTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchSymmTensorField,
    processorUnallocatedFvsPatchSymmTensorField,
    patch
);

// Tensor
defineNamedTemplateTypeNameAndDebug(processorUnallocatedFvsPatchTensorField, 0);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchTensorField,
    processorUnallocatedFvsPatchTensorField,
    dictionary
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchTensorField,
    processorUnallocatedFvsPatchTensorField,
    patchMapper
);
addToRunTimeSelectionTable
(
    unallocatedFvsPatchTensorField,
    processorUnallocatedFvsPatchTensorField,
    patch
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
