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

\*---------------------------------------------------------------------------*/

#include "unallocatedFvMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(unallocatedFvMesh, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unallocatedFvMesh::unallocatedFvMesh
(
    //const fvMesh& procMesh,
    const objectRegistry& db,
    const label nCells,
    const unallocatedFvBoundaryMesh& boundary,
    const globalMeshData& globalData
)
:
    volMesh(*this), //(procMesh),
    db_(db),
    nCells_(nCells),
    boundary_(boundary),
    globalData_(globalData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::unallocatedFvMesh::~unallocatedFvMesh()
{
    if (debug)
    {
        Pout<< "~unallocatedFvMesh::unallocatedFvMesh()"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::unallocatedFvMesh::operator!=(const unallocatedFvMesh& bm) const
{
    return &bm != this;
}


bool Foam::unallocatedFvMesh::operator==(const unallocatedFvMesh& bm) const
{
    return &bm == this;
}


// ************************************************************************* //
