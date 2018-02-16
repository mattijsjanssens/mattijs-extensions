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

#include "uVolFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplate2TypeNameAndDebug(uVolScalarField::Internal, 0);
// defineTemplate2TypeNameAndDebug(uVolVectorField::Internal, 0);
// defineTemplate2TypeNameAndDebug
// (
//     uVolSphericalTensorField::Internal,
//     0
// );
// defineTemplate2TypeNameAndDebug
// (
//     uVolSymmTensorField::Internal,
//     0
// );
// defineTemplate2TypeNameAndDebug(uVolTensorField::Internal, 0);

defineTemplateTypeNameAndDebugWithName(uVolScalarField, "volScalarField", 0);
// defineTemplateTypeNameAndDebugWithName(uVolVectorField, "volVectorField", 0);
// defineTemplateTypeNameAndDebugWithName
// (
//     uVolSphericalTensorField,
//     "volSphericalTensorField",
//     0
// );
// defineTemplateTypeNameAndDebugWithName
// (
//     uVolSymmTensorField,
//     "volSymmTensorField",
//     0
// );
// defineTemplateTypeNameAndDebugWithName(uVolTensorField, "volTensorField", 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// specialization for scalar fields
template<>
tmp<GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>>
GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>::component
(
    const direction
) const
{
    return *this;
}


// specialization for scalar fields
template<>
void GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>::replace
(
    const direction,
    const GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>&
    gsf
)
{
    *this == gsf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
