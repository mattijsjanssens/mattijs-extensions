/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "scalarField.H"
#include "GAMGInterface.H"
#include "GAMGAgglomeration.H"
#include "lduMatrix.H"


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::GAMGInterface> Foam::GAMGInterface::New
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
{
    const word coupleType(fineInterface.type());

    auto cstrIter = lduInterfaceConstructorTablePtr_->cfind(coupleType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown GAMGInterface type " << coupleType << ".\n"
            << "Valid GAMGInterface types :"
            << lduInterfaceConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<GAMGInterface>
    (
        cstrIter()
        (
            index,
            coarseInterfaces,
            fineInterface,
            localRestrictAddressing,
            neighbourRestrictAddressing,
            fineLevelIndex,
            coarseComm
        )
    );
}


Foam::autoPtr<Foam::GAMGInterface> Foam::GAMGInterface::New
(
    const word& coupleType,
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
{
    auto cstrIter = IstreamConstructorTablePtr_->cfind(coupleType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown GAMGInterface type " << coupleType << ".\n"
            << "Valid GAMGInterface types :"
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<GAMGInterface>(cstrIter()(index, coarseInterfaces, is));
}


// ************************************************************************* //
