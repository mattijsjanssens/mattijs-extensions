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

#include "surfaceFieldStreamReconstructor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// No need for typedefs outside this usage so no need for .H file
typedef surfaceFieldStreamReconstructor<scalar>
surfaceScalarFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    streamReconstructor,
    surfaceScalarFieldStreamReconstructor,
    word,
    surfaceScalarField
);

typedef surfaceFieldStreamReconstructor<vector>
surfaceVectorFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    streamReconstructor,
    surfaceVectorFieldStreamReconstructor,
    word,
    surfaceVectorField
);

typedef surfaceFieldStreamReconstructor<sphericalTensor>
surfaceSphericalTensorFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    streamReconstructor,
    surfaceSphericalTensorFieldStreamReconstructor,
    word,
    surfaceSphericalTensorField
);

typedef surfaceFieldStreamReconstructor<symmTensor>
surfaceSymmTensorFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    streamReconstructor,
    surfaceSymmTensorFieldStreamReconstructor,
    word,
    surfaceSymmTensorField
);

typedef surfaceFieldStreamReconstructor<tensor>
surfaceTensorFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    streamReconstructor,
    surfaceTensorFieldStreamReconstructor,
    word,
    surfaceTensorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
