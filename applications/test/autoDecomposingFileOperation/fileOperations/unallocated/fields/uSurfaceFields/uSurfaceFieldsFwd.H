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

InClass
    Foam::uSurfaceFields

Description

SourceFiles
    uSurfaceFields.C

\*---------------------------------------------------------------------------*/

#ifndef uSurfaceFieldsFwd_H
#define uSurfaceFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class unallocatedSurfaceMesh;

template<class Type>
class unallocatedFvsPatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

typedef GeometricField<scalar, unallocatedFvsPatchField, unallocatedSurfaceMesh>
uSurfaceScalarField;
typedef GeometricField<vector, unallocatedFvsPatchField, unallocatedSurfaceMesh>
uSurfaceVectorField;
typedef
GeometricField
<
    sphericalTensor,
    unallocatedFvsPatchField,
    unallocatedSurfaceMesh
>
uSurfaceSphericalTensorField;
typedef GeometricField
<
    symmTensor,
    unallocatedFvsPatchField,
    unallocatedSurfaceMesh
>
uSurfaceSymmTensorField;
typedef GeometricField<tensor, unallocatedFvsPatchField, unallocatedSurfaceMesh>
uSurfaceTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
