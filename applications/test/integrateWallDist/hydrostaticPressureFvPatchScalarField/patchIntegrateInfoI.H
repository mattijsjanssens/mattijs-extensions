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

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::patchIntegrateInfo::patchIntegrateInfo()
:
    value_(-GREAT)
{}


// Construct from origin, distance
inline Foam::patchIntegrateInfo::patchIntegrateInfo(const scalar value)
:
    value_(value)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::patchIntegrateInfo::value() const
{
    return value_;
}


template<class TrackingData>
inline bool Foam::patchIntegrateInfo::valid(TrackingData& td) const
{
    return value_ != -GREAT;
}


template<class TrackingData>
inline void Foam::patchIntegrateInfo::transform
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const tensor& rotTensor,
    const scalar tol,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::patchIntegrateInfo::updateEdge
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const label edgei,
    const label facei,
    const patchIntegrateInfo& faceInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Update edge from face value
    if (!valid(td))
    {
        scalar dh = td.edgeCentres_[edgei]-td.faceCentres_[facei];
        scalar drho = 0.5*(td.fFld_[facei]+td.eFld_[edgei]);
        value_ = faceInfo.value_ + dh*drho;
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::patchIntegrateInfo::updateEdge
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const patchIntegrateInfo& edgeInfo,
    const bool sameOrientation,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        value_ = edgeInfo.value_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::patchIntegrateInfo::updateFace
(
    const polyMesh& mesh,
    const primitivePatch& patch,
    const label facei,
    const label edgei,
    const patchIntegrateInfo& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    // Update face from edge value
    if (!valid(td))
    {
        scalar dh = td.faceCentres_[facei]-td.edgeCentres_[edgei];
        scalar drho = 0.5*(td.fFld_[facei]+td.eFld_[edgei]);
        value_ = edgeInfo.value_ + dh*drho;
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::patchIntegrateInfo::equal
(
    const patchIntegrateInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::patchIntegrateInfo::operator==
(
    const Foam::patchIntegrateInfo& rhs
) const
{
    return value() == rhs.value();
}


inline bool Foam::patchIntegrateInfo::operator!=
(
    const Foam::patchIntegrateInfo& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::patchIntegrateInfo& wDist
)
{
    return os << wDist.value();
}


inline Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::patchIntegrateInfo& wDist
)
{
    return is >> wDist.value_;
}


// ************************************************************************* //
