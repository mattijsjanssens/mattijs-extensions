/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "lduPrimitiveProcessorInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduPrimitiveProcessorInterface, 0);
    addToRunTimeSelectionTable
    (
        lduInterface,
        lduPrimitiveProcessorInterface,
        Istream
    );
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveProcessorInterface::lduPrimitiveProcessorInterface
(
    const labelUList& faceCells,
    const label comm,
    const label myProcNo,
    const label neighbProcNo,
    const tensorField& forwardT,
    const int tag
)
:
    faceCells_(faceCells),
    comm_(comm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    forwardT_(forwardT),
    tag_(tag)
{}


Foam::lduPrimitiveProcessorInterface::lduPrimitiveProcessorInterface
(
    Istream& is
)
:
    faceCells_(is),
    comm_(readLabel(is)),
    myProcNo_(readLabel(is)),
    neighbProcNo_(readLabel(is)),
    forwardT_(is),
    tag_(readInt(is))
{}


Foam::lduPrimitiveProcessorInterface::lduPrimitiveProcessorInterface
(
    const lduPrimitiveProcessorInterface& pp,
    const label index,
    const labelUList& mapAddressing
)
:
    faceCells_(UIndirectList<label>(pp.faceCells_, mapAddressing)),
    comm_(pp.comm_),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    forwardT_(pp.forwardT_),
    tag_(pp.tag_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lduPrimitiveProcessorInterface::~lduPrimitiveProcessorInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lduPrimitiveProcessorInterface::write(Ostream& os) const
{
    os  << faceCells_ << token::SPACE
        << comm_ << token::SPACE
        << myProcNo_ << token::SPACE
        << neighbProcNo_ << token::SPACE
        << forwardT_ << token::SPACE
        << tag_;
}


// ************************************************************************* //
