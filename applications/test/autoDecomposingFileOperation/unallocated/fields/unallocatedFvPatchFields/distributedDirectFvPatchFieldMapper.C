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

#include "distributedDirectFvPatchFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedDirectFvPatchFieldMapper::distributedDirectFvPatchFieldMapper
(
    const labelUList& directAddressing,
    const mapDistributeBase& distMap,
    const label constructSize
)
:
    directAddressing_(directAddressing),
    hasUnmapped_(false)
{
    if (constructSize == -1)
    {
        distMap_ = distMap;
    }
    else
    {
        distMap_.constructSize() = constructSize;
        distMap_.subMap() = distMap.constructMap();
        distMap_.subHasFlip() = distMap.constructHasFlip();
        distMap_.constructMap() = distMap.subMap();
        distMap_.constructHasFlip() = distMap.subHasFlip();
    }

DebugVar(distMap_.constructSize());
DebugVar(distMap_.subMap());
DebugVar(distMap_.constructMap());

    if
    (
        notNull(directAddressing_)
     && directAddressing_.size()
     && min(directAddressing_) < 0
    )
    {
        hasUnmapped_ = true;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
