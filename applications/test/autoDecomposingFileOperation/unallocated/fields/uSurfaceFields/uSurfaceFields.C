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

#include "uSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplate2TypeNameAndDebug(uSurfaceScalarField::Internal, 0);
defineTemplate2TypeNameAndDebug(uSurfaceVectorField::Internal, 0);
defineTemplate2TypeNameAndDebug(uSurfaceSphericalTensorField::Internal, 0);
defineTemplate2TypeNameAndDebug(uSurfaceSymmTensorField::Internal, 0);
defineTemplate2TypeNameAndDebug(uSurfaceTensorField::Internal, 0);

defineTemplateTypeNameAndDebugWithName
(
    uSurfaceScalarField,
    "surfaceScalarField",
    0
);
defineTemplateTypeNameAndDebugWithName
(
    uSurfaceVectorField,
    "surfaceVectorField",
    0
);
defineTemplateTypeNameAndDebugWithName
(
    uSurfaceSphericalTensorField,
    "surfaceSphericalTensorField",
    0
);
defineTemplateTypeNameAndDebugWithName
(
    uSurfaceSymmTensorField,
    "surfaceSymmTensorField",
    0
);
defineTemplateTypeNameAndDebugWithName
(
    uSurfaceTensorField,
    "surfaceTensorField",
    0
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// specialization for scalar fields
template<>
tmp<GeometricField<scalar, unallocatedFvsPatchField, unallocatedSurfaceMesh>>
GeometricField<scalar, unallocatedFvsPatchField, unallocatedSurfaceMesh>::
component
(
    const direction
) const
{
    return *this;
}


// specialization for scalar fields
template<>
void GeometricField<scalar, unallocatedFvsPatchField, unallocatedSurfaceMesh>::
replace
(
    const direction,
    const GeometricField
    <
        scalar,
        unallocatedFvsPatchField,
        unallocatedSurfaceMesh
    >&
    gsf
)
{
    *this == gsf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
