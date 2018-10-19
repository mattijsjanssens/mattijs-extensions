/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(unallocatedFvMesh, 0);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unallocatedFvMesh::unallocatedFvMesh
(
    const word& name,
    const objectRegistry& db,
    const label nCells,
    const wordList& patchNames,
    const wordList& patchTypes,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const labelListList& basePatchAddressing,
    const List<wordList>& basePatchGroups,
    const PtrList<dictionary>& patchDicts,
    const globalMeshData& globalData
)
:
    name_(name),
    db_(db),
    nInternalFaces_(patchStarts.size() ? min(patchStarts) : 0),
    nCells_(nCells),
    globalData_(globalData)
{
    boundary_.setSize(patchNames.size());
    forAll(boundary_, patchi)
    {
        boundary_.set
        (
            patchi,
            new unallocatedGenericFvPatch
            (
                *reinterpret_cast<const polyPatch*>(0), //unused reference
                patchNames[patchi],
                patchTypes[patchi],
                patchSizes[patchi],
                patchStarts[patchi],
                basePatchAddressing[patchi],
                patchDicts[patchi],
                patchi,
                *reinterpret_cast<const fvBoundaryMesh*>(0),// unused reference
                basePatchGroups[patchi]
            )
        );
    }

    nFaces_ = 0;
    forAll(patchStarts, patchi)
    {
        nFaces_ = max(nFaces_, patchStarts[patchi]+patchSizes[patchi]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::unallocatedFvMesh::~unallocatedFvMesh()
{
    if (debug)
    {
        Pout<< indent
            << "~unallocatedFvMesh::unallocatedFvMesh()" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::unallocatedFvMesh::operator!=(const unallocatedFvMesh& bm) const
{
    return &bm != this;
}


bool Foam::unallocatedFvMesh::operator==(const unallocatedFvMesh& bm) const
{
    return &bm == this;
}


// * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<unallocatedFvMesh>& proxy
)
{
    const unallocatedFvMesh& mesh = proxy.t_;

    const unallocatedFvBoundaryMesh& patches = mesh.boundary();

    os.setf(ios_base::left);

    os  << indent << setw(20) << "name" << mesh.name() << nl
        << indent << setw(20) << "faces" << mesh.nFaces() << nl
        << indent << setw(20) << "internal faces" << mesh.nInternalFaces() << nl
        << indent << setw(20) << "cells" << mesh.nCells() << nl
        << indent << setw(20) << "boundary patches" << patches.size() << nl
        << nl;

    os.setf(ios_base::left);

    os  << indent
        << setw(20) << "Patch"
        << setw(9) << "Faces"
        << setw(9) << "Type"
        << nl;
    forAll(patches, patchi)
    {
        const unallocatedGenericFvPatch& pp = patches[patchi];
        os  << indent
            << setw(20) << pp.name()
            << setw(9) << pp.size()
            << setw(9) << pp.type()
            << nl;
    }

    return os;
}


// ************************************************************************* //
