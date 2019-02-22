/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "lduPrimitiveInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduPrimitiveInterface, 0);
    addToRunTimeSelectionTable(lduInterface, lduPrimitiveInterface, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveInterface::lduPrimitiveInterface
(
    const labelUList& faceCells,
    const labelUList& faceGlobalCells,
    const labelUList& nbrFaceGlobalCells
)
:
    faceCells_(faceCells),
    faceGlobalCells_(faceGlobalCells),
    nbrFaceGlobalCells_(nbrFaceGlobalCells)
{}


Foam::lduPrimitiveInterface::lduPrimitiveInterface(Istream& is)
:
    faceCells_(is),
    faceGlobalCells_(is),
    nbrFaceGlobalCells_(is)
{}


Foam::lduPrimitiveInterface::lduPrimitiveInterface
(
    const lduPrimitiveInterface& pp,
    const label index,
    const labelUList& mapAddressing
)
:
    faceCells_(UIndirectList<label>(pp.faceCells_, mapAddressing)),
    faceGlobalCells_(UIndirectList<label>(pp.faceGlobalCells_, mapAddressing)),
    nbrFaceGlobalCells_
    (
        UIndirectList<label>(pp.nbrFaceGlobalCells_, mapAddressing)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lduPrimitiveInterface::~lduPrimitiveInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lduPrimitiveInterface::write(Ostream& os) const
{
    os  << faceCells_ << token::SPACE
        << faceGlobalCells_ << token::SPACE
        << nbrFaceGlobalCells_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
