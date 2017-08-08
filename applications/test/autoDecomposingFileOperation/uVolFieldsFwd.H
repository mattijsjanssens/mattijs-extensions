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

InClass
    Foam::uVolFields

Description

SourceFiles
    uVolFields.C

\*---------------------------------------------------------------------------*/

#ifndef uVolFieldsFwd_H
#define uVolFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class unallocatedFvMesh;

template<class Type>
class unallocatedFvPatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

typedef GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>
uVolScalarField;
typedef GeometricField<vector, unallocatedFvPatchField, unallocatedFvMesh>
uVolVectorField;
typedef
GeometricField<sphericalTensor, unallocatedFvPatchField, unallocatedFvMesh>
uVolSphericalTensorField;
typedef GeometricField<symmTensor, unallocatedFvPatchField, unallocatedFvMesh>
uVolSymmTensorField;
typedef GeometricField<tensor, unallocatedFvPatchField, unallocatedFvMesh>
uVolTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
