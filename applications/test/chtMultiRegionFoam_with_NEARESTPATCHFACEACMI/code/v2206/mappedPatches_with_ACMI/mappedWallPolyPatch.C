/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "mappedWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedWallPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mappedWallPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedWallPolyPatch::mappedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, size, start, index, bm, patchType),
    mappedPatchBase2(static_cast<const polyPatch&>(*this))
{
    //  mapped is not constraint type so add mapped group explicitly
    inGroups().appendUniq(mappedPolyPatch2::typeName);
}


Foam::mappedWallPolyPatch::mappedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const mappedPatchBase2::sampleMode mode,
    const word& samplePatch,
    const vectorField& offset,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(name, size, start, index, bm, typeName),
    mappedPatchBase2
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::mappedWallPolyPatch::mappedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const mappedPatchBase2::sampleMode mode,
    const word& samplePatch,
    const vector& offset,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(name, size, start, index, bm, typeName),
    mappedPatchBase2
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::mappedWallPolyPatch::mappedWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, dict, index, bm, patchType),
    mappedPatchBase2(*this, dict)
{
    //  mapped is not constraint type so add mapped group explicitly
    inGroups().appendUniq(mappedPolyPatch2::typeName);
}


Foam::mappedWallPolyPatch::mappedWallPolyPatch
(
    const mappedWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(pp, bm),
    mappedPatchBase2(*this, pp)
{}


Foam::mappedWallPolyPatch::mappedWallPolyPatch
(
    const mappedWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, newSize, newStart),
    mappedPatchBase2(*this, pp)
{}


Foam::mappedWallPolyPatch::mappedWallPolyPatch
(
    const mappedWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, mapAddressing, newStart),
    mappedPatchBase2(*this, pp, mapAddressing)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedWallPolyPatch::~mappedWallPolyPatch()
{
    mappedPatchBase2::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedWallPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::initGeometry(pBufs);
}


void Foam::mappedWallPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    wallPolyPatch::calcGeometry(pBufs);
    mappedPatchBase2::clearOut();
}


void Foam::mappedWallPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::initMovePoints(pBufs, p);
}


void Foam::mappedWallPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    wallPolyPatch::movePoints(pBufs, p);
    mappedPatchBase2::clearOut();
}


void Foam::mappedWallPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    wallPolyPatch::initUpdateMesh(pBufs);
}


void Foam::mappedWallPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    wallPolyPatch::updateMesh(pBufs);
    mappedPatchBase2::clearOut();
}


void Foam::mappedWallPolyPatch::write(Ostream& os) const
{
    wallPolyPatch::write(os);
    mappedPatchBase2::write(os);
}


// ************************************************************************* //
