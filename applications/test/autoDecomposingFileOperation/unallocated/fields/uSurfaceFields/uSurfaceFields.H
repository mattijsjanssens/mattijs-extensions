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

#ifndef uSurfaceFields_H
#define uSurfaceFields_H

#include "GeometricFields.H"
#include "unallocatedSurfaceMesh.H"
#include "unallocatedFvsPatchField.H"
//#include "fvMesh.H"
//#include "fvsPatchField.H"
#include "uSurfaceFieldsFwd.H"
//#include "calculatedFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<GeometricField<scalar, unallocatedFvsPatchField, unallocatedSurfaceMesh>>
GeometricField<scalar, unallocatedFvsPatchField, unallocatedSurfaceMesh>::
component
(
    const direction
) const;

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
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "uSurfaceFieldsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
