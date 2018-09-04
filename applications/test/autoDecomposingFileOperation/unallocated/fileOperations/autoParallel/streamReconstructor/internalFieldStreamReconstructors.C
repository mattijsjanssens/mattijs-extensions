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

#include "internalFieldStreamReconstructor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//- Version of addToRunTimeSelectionTable with explicitly provided
//  object name such that we can use e.g. 'volScalarField::Internal'
//  as a lookup
#define addExplicitNamedToRunTimeSelectionTable\
(baseType,thisType,argNames,addName,lookup)                                    \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType::add##argNames##ConstructorToTable<thisType>                      \
        addName(#lookup)


// No need for typedefs outside this usage so no need for .H file

typedef internalFieldStreamReconstructor<scalar>
internalScalarFieldStreamReconstructor;
addExplicitNamedToRunTimeSelectionTable
(
    streamReconstructor,
    internalScalarFieldStreamReconstructor,
    word,
    add_volScalarFieldInternal,
    volScalarField::Internal
);

typedef internalFieldStreamReconstructor<vector>
internalVectorFieldStreamReconstructor;
addExplicitNamedToRunTimeSelectionTable
(
    streamReconstructor,
    internalVectorFieldStreamReconstructor,
    word,
    add_volVectorFieldInternal,
    volVectorField::Internal
);

typedef internalFieldStreamReconstructor<sphericalTensor>
internalSphericalTensorFieldStreamReconstructor;
addExplicitNamedToRunTimeSelectionTable
(
    streamReconstructor,
    internalSphericalTensorFieldStreamReconstructor,
    word,
    add_volSphericalTensorFieldInternal,
    volSphericalTensorField::Internal
);

typedef internalFieldStreamReconstructor<symmTensor>
internalSymmTensorFieldStreamReconstructor;
addExplicitNamedToRunTimeSelectionTable
(
    streamReconstructor,
    internalSymmTensorFieldStreamReconstructor,
    word,
    add_volSymmTensorFieldInternal,
    volSymmTensorField::Internal
);

typedef internalFieldStreamReconstructor<tensor>
internalTensorFieldStreamReconstructor;
addExplicitNamedToRunTimeSelectionTable
(
    streamReconstructor,
    internalTensorFieldStreamReconstructor,
    word,
    add_volTensorFieldInternal,
    volTensorField::Internal
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
