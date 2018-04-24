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

#include "unallocatedFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(unallocatedFvsPatchScalarField, 0);
defineTemplateRunTimeSelectionTable(unallocatedFvsPatchScalarField, patch);
defineTemplateRunTimeSelectionTable(unallocatedFvsPatchScalarField, dictionary);
defineTemplateRunTimeSelectionTable
(
    unallocatedFvsPatchScalarField,
    patchMapper
);

defineNamedTemplateTypeNameAndDebug(unallocatedFvsPatchVectorField, 0);
defineTemplateRunTimeSelectionTable(unallocatedFvsPatchVectorField, patch);
defineTemplateRunTimeSelectionTable(unallocatedFvsPatchVectorField, dictionary);
defineTemplateRunTimeSelectionTable
(
    unallocatedFvsPatchVectorField,
    patchMapper
);

defineNamedTemplateTypeNameAndDebug
(
    unallocatedFvsPatchSphericalTensorField,
    0
);
defineTemplateRunTimeSelectionTable
(
    unallocatedFvsPatchSphericalTensorField,
    patch
);
defineTemplateRunTimeSelectionTable
(
    unallocatedFvsPatchSphericalTensorField,
    dictionary
);
defineTemplateRunTimeSelectionTable
(
    unallocatedFvsPatchSphericalTensorField,
    patchMapper
);

defineNamedTemplateTypeNameAndDebug(unallocatedFvsPatchSymmTensorField, 0);
defineTemplateRunTimeSelectionTable(unallocatedFvsPatchSymmTensorField, patch);
defineTemplateRunTimeSelectionTable
(
    unallocatedFvsPatchSymmTensorField,
    dictionary
);
defineTemplateRunTimeSelectionTable
(
    unallocatedFvsPatchSymmTensorField,
    patchMapper
);

defineNamedTemplateTypeNameAndDebug(unallocatedFvsPatchTensorField, 0);
defineTemplateRunTimeSelectionTable(unallocatedFvsPatchTensorField, patch);
defineTemplateRunTimeSelectionTable(unallocatedFvsPatchTensorField, dictionary);
defineTemplateRunTimeSelectionTable
(
    unallocatedFvsPatchTensorField,
    patchMapper
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
