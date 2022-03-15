/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 M. Janssens
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

template<class Type>
inline Foam::limitedDistanceData<Type>::limitedDistanceData()
{}


template<class Type>
inline Foam::limitedDistanceData<Type>::limitedDistanceData
(
    const List<scalar>& weights,
    const List<scalar>& maxDist,
    const List<point>& origin,
    const List<Type>& data
)
:
    weights_(weights),
    maxDist_(maxDist),
    origin_(origin),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::scalar Foam::limitedDistanceData<Type>::weightFunction::weight
(
    const polyMesh& mesh,
    const label celli,
    const point& origin,
    const scalar r
) const
{
    // Rough approximation: count number of points inside sphere

    const scalar d(mag(mesh.cellCentres()[celli]-origin));
    const scalar scale(fun_.value(d/r));
    if (scale < ROOTVSMALL)
    {
        return 0.0;
    }
    else
    {
        label nInside = 0;

        const labelList& cPoints = mesh.cellPoints(celli);
        for (const label pointi : cPoints)
        {
            const point& pt = mesh.points()[pointi];
            //if (magSqr(pt-origin) < r2)
            if (fun_.value(mag(pt-origin)/r) >= ROOTVSMALL)
            {
                nInside++;
            }
        }

        return scale*mesh.cellVolumes()[celli]*nInside/cPoints.size();
    }
}


template<class Type>
template<class TrackingData>
inline bool Foam::limitedDistanceData<Type>::valid(TrackingData& td) const
{
//    return distance_ != -1;
    return true;
}


// No geometric data so never any problem on cyclics
template<class Type>
template<class TrackingData>
inline bool Foam::limitedDistanceData<Type>::sameGeometry
(
    const polyMesh&,
    const limitedDistanceData<Type>&,
    const scalar,
    TrackingData&
) const
{
    return true;
}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::limitedDistanceData<Type>::leaveDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData&
)
{
    for (auto& o : origin_)
    {
        o -= faceCentre;
    }
}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::limitedDistanceData<Type>::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData&
)
{
    for (auto& o : origin_)
    {
        o = Foam::transform(rotTensor, o);
    }
}


// No geometric data.
template<class Type>
template<class TrackingData>
inline void Foam::limitedDistanceData<Type>::enterDomain
(
    const polyMesh&,
    const polyPatch& patch,
    const label patchFacei,
    const point& faceCentre,
    TrackingData&
)
{
    // back to absolute form
    for (auto& o : origin_)
    {
        o += faceCentre;
    }
}


// Update cell with neighbouring face information
template<class Type>
template<class TrackingData>
inline bool Foam::limitedDistanceData<Type>::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const limitedDistanceData<Type>& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    bool changed = false;

    forAll(neighbourInfo.data_, i)
    {
        const label index = data_.find(neighbourInfo.data_[i]);
        if (index == -1)
        {
            const scalar d = neighbourInfo.maxDist_[i];
            const point& o = neighbourInfo.origin_[i];

            // Check if cell is inside search sphere
            //const point& pt = mesh.cellCentres()[thisCelli];
            // Check if face is inside search sphere
            //const point& pt = mesh.faceCentres()[neighbourFacei];
            //if (magSqr(pt-o) < d)

            const scalar weight = td.weight(mesh, thisCelli, o, d);

            //Pout<< "For cell:" << thisCelli
            //    << " at:" << mesh.cellCentres()[thisCelli]
            //    << " have particle origin:" << o
            //    << ", d:" << d
            //    << " overlap:" << overlapVol
            //    << " coming from face:"
            //    << mesh.faceCentres()[neighbourFacei]
            //    << endl;



            if (weight > ROOTVSMALL)
            {

                maxDist_.append(d);
                origin_.append(o);
                data_.append(neighbourInfo.data_[i]);
                weights_.append(weight);
                changed = true;
            }
        }
    }

    return changed;
}


// Update face with neighbouring cell information
template<class Type>
template<class TrackingData>
inline bool Foam::limitedDistanceData<Type>::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const limitedDistanceData<Type>& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    // From cell to its faces.

    bool changed = false;

    forAll(neighbourInfo.data_, i)
    {
        const label index = data_.find(neighbourInfo.data_[i]);
        if (index == -1)
        {
            //const scalar nbrD = neighbourInfo.maxDist_[i];
            // Check if cell is inside search sphere
            //const point& pt = mesh.cellCentres()[neighbourCelli];
            // Check if face is inside search sphere
            //const point& pt = mesh.faceCentres()[thisFacei];
            //if (magSqr(pt-neighbourInfo.origin_[i]) < nbrD)
            {
                maxDist_.append(neighbourInfo.maxDist_[i]);
                origin_.append(neighbourInfo.origin_[i]);
                data_.append(neighbourInfo.data_[i]);

                // Calculate overlap volume?
                weights_.append(0);

                changed = true;
            }
        }
    }

    return changed;
}


// Update face with coupled face information
template<class Type>
template<class TrackingData>
inline bool Foam::limitedDistanceData<Type>::updateFace
(
    const polyMesh&,
    const label thisFacei,
    const limitedDistanceData<Type>& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    // From face to face (e.g. coupled faces)
    bool changed = false;

    forAll(neighbourInfo.data_, i)
    {
        const label index = data_.find(neighbourInfo.data_[i]);
        if (index == -1)
        {
            maxDist_.append(neighbourInfo.maxDist_[i]);
            origin_.append(neighbourInfo.origin_[i]);
            data_.append(neighbourInfo.data_[i]);

            // Calculate overlap volume?
            weights_.append(0);

            changed = true;
        }
    }

    return changed;
}


template<class Type>
template<class TrackingData>
inline bool Foam::limitedDistanceData<Type>::equal
(
    const limitedDistanceData<Type>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::limitedDistanceData<Type>::operator==
(
    const Foam::limitedDistanceData<Type>& rhs
) const
{
    return data_ == rhs.data_;
}


template<class Type>
inline bool Foam::limitedDistanceData<Type>::operator!=
(
    const Foam::limitedDistanceData<Type>& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
