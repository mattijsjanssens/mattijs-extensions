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

#include "ACMIWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedPolyPatch.H"
#include "cyclicACMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ACMIWallPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, ACMIWallPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        ACMIWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::ACMIWallPolyPatch::ACMIWallPolyPatch
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
    coupleGroup_(),
    ACMIPatch_(""),
    ACMIPatchID_(-1)
{}


Foam::ACMIWallPolyPatch::ACMIWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    wallPolyPatch(name, dict, index, bm, patchType),
    coupleGroup_(dict),
    ACMIPatch_(dict.lookupOrDefault<word>("ACMIPatch", "")),
    ACMIPatchID_(-1)
{}


Foam::ACMIWallPolyPatch::ACMIWallPolyPatch
(
    const ACMIWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(pp, bm),
    coupleGroup_(pp.coupleGroup_),
    ACMIPatch_(pp.ACMIPatch_),
    ACMIPatchID_(-1)
{}


Foam::ACMIWallPolyPatch::ACMIWallPolyPatch
(
    const ACMIWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, newSize, newStart),
    coupleGroup_(pp.coupleGroup_),
    ACMIPatch_(pp.ACMIPatch_),
    ACMIPatchID_(-1)
{}


Foam::ACMIWallPolyPatch::ACMIWallPolyPatch
(
    const ACMIWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, mapAddressing, newStart),
    coupleGroup_(pp.coupleGroup_),
    ACMIPatch_(pp.ACMIPatch_),
    ACMIPatchID_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ACMIWallPolyPatch::~ACMIWallPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//void Foam::ACMIWallPolyPatch::initGeometry(PstreamBuffers& pBufs)
//{
//    wallPolyPatch::initGeometry(pBufs);
//}
//
//
//void Foam::ACMIWallPolyPatch::calcGeometry(PstreamBuffers& pBufs)
//{
//    wallPolyPatch::calcGeometry(pBufs);
//}
//
//
//void Foam::ACMIWallPolyPatch::initMovePoints
//(
//    PstreamBuffers& pBufs,
//    const pointField& p
//)
//{
//    wallPolyPatch::initMovePoints(pBufs, p);
//}
//
//
//void Foam::ACMIWallPolyPatch::movePoints
//(
//    PstreamBuffers& pBufs,
//    const pointField& p
//)
//{
//    wallPolyPatch::movePoints(pBufs, p);
//}
//
//
//void Foam::ACMIWallPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
//{
//    wallPolyPatch::initUpdateMesh(pBufs);
//}
//
//
//void Foam::ACMIWallPolyPatch::updateMesh(PstreamBuffers& pBufs)
//{
//    wallPolyPatch::updateMesh(pBufs);
//}


const Foam::cyclicACMIPolyPatch& Foam::ACMIWallPolyPatch::ACMI() const
{
    if (ACMIPatchID_ == -1)
    {
        if (!ACMIPatch_.empty())
        {
            ACMIPatchID_ = boundaryMesh().findPatchID(ACMIPatch_);
        }
        else
        {
            // Use couple group
            ACMIPatchID_ = coupleGroup_.findOtherPatchID(*this);
        }
    }
    return refCast<const cyclicACMIPolyPatch>(boundaryMesh()[ACMIPatchID_]);
}


void Foam::ACMIWallPolyPatch::write(Ostream& os) const
{
    wallPolyPatch::write(os);
    coupleGroup_.write(os);
    if (!ACMIPatch_.empty())
    {
        os.writeKeyword("ACMIPatch") << ACMIPatch_
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
