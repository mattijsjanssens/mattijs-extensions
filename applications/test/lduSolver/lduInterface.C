/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "lduInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduInterface, 0);
    defineRunTimeSelectionTable(lduInterface, dictionary);
    defineRunTimeSelectionTable(lduInterface, Istream);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lduInterface::~lduInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::lduInterface> Foam::lduInterface::New
(
    const dictionary& dict,
    const label index
)
{
    word patchType(dict.lookup("type"));

    if (debug)
    {
        InfoInFunction << "Constructing lduInterface " << patchType << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(patchType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown lduInterface type "
            << patchType << nl << nl
            << "Valid lduInterface types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<lduInterface>(cstrIter()(dict, index));
}


Foam::autoPtr<Foam::lduInterface> Foam::lduInterface::New
(
    const word& patchType,
    Istream& is
)
{
    if (debug)
    {
        InfoInFunction << "Constructing lduInterface " << patchType << endl;
    }

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(patchType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Unknown lduInterface type "
            << patchType << nl << nl
            << "Valid lduInterface types are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<lduInterface>(cstrIter()(is));
}


// ************************************************************************* //
