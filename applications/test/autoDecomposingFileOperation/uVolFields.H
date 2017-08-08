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

#ifndef uVolFields_H
#define uVolFields_H

#include "GeometricFields.H"
#include "unallocatedFvMesh.H"
#include "unallocatedFvPatchField.H"
//#include "fvMesh.H"
//#include "fvPatchField.H"
#include "uVolFieldsFwd.H"
//#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>>
GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>::component
(
    const direction
) const;

template<>
void GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>::
replace
(
    const direction,
    const GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>& sf
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "uVolFieldsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
