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

#include "FieldStreamReconstructor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// No need for typedefs outside this usage so no need for .H file
typedef FieldStreamReconstructor<label> labelFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    lagrangianStreamReconstructor,
    labelFieldStreamReconstructor,
    cloudName,
    labelField
);

typedef FieldStreamReconstructor<scalar> scalarFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    lagrangianStreamReconstructor,
    scalarFieldStreamReconstructor,
    cloudName,
    scalarField
);

typedef FieldStreamReconstructor<vector> vectorFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    lagrangianStreamReconstructor,
    vectorFieldStreamReconstructor,
    cloudName,
    vectorField
);

typedef FieldStreamReconstructor<sphericalTensor>
sphericalTensorFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    lagrangianStreamReconstructor,
    sphericalTensorFieldStreamReconstructor,
    cloudName,
    sphericalTensorField
);

typedef FieldStreamReconstructor<symmTensor>
symmTensorFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    lagrangianStreamReconstructor,
    symmTensorFieldStreamReconstructor,
    cloudName,
    symmTensorField
);

typedef FieldStreamReconstructor<tensor> tensorFieldStreamReconstructor;
addNamedToRunTimeSelectionTable
(
    lagrangianStreamReconstructor,
    tensorFieldStreamReconstructor,
    cloudName,
    tensorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
