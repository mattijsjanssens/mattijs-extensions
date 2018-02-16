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

#include "unallocatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(unallocatedFvPatchScalarField, 0);
defineTemplateRunTimeSelectionTable(unallocatedFvPatchScalarField, patch);
defineTemplateRunTimeSelectionTable(unallocatedFvPatchScalarField, dictionary);
defineTemplateRunTimeSelectionTable(unallocatedFvPatchScalarField, patchMapper);

// defineNamedTemplateTypeNameAndDebug(unallocatedFvPatchVectorField, 0);
// defineTemplateRunTimeSelectionTable(unallocatedFvPatchVectorField, patch);
// defineTemplateRunTimeSelectionTable(unallocatedFvPatchVectorField, dictionary);
// 
// defineNamedTemplateTypeNameAndDebug
// (
//     unallocatedFvPatchSphericalTensorField,
//     0
// );
// defineTemplateRunTimeSelectionTable
// (
//     unallocatedFvPatchSphericalTensorField,
//     patch
// );
// defineTemplateRunTimeSelectionTable
// (
//     unallocatedFvPatchSphericalTensorField,
//     dictionary
// );
// 
// defineNamedTemplateTypeNameAndDebug(unallocatedFvPatchSymmTensorField, 0);
// defineTemplateRunTimeSelectionTable(unallocatedFvPatchSymmTensorField, patch);
// defineTemplateRunTimeSelectionTable
// (
//     unallocatedFvPatchSymmTensorField,
//     dictionary
// );
// 
// defineNamedTemplateTypeNameAndDebug(unallocatedFvPatchTensorField, 0);
// defineTemplateRunTimeSelectionTable(unallocatedFvPatchTensorField, patch);
// defineTemplateRunTimeSelectionTable(unallocatedFvPatchTensorField, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
