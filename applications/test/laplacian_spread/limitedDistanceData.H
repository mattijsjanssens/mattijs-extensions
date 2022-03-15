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

Class
    Foam::limitedDistanceData

Description
    For use with FaceCellWave. Determines topological distance to starting
    faces. Templated on passive transported data.

SourceFiles
    limitedDistanceDataI.H
    limitedDistanceData.C

\*---------------------------------------------------------------------------*/

#ifndef limitedDistanceData_H
#define limitedDistanceData_H

#include "point.H"
#include "tensor.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;
class polyMesh;
template<class Type> class limitedDistanceData;

template<class Type>
Istream& operator>>(Istream&, limitedDistanceData<Type>&);
template<class Type>
Ostream& operator<<(Ostream&, const limitedDistanceData<Type>&);

/*---------------------------------------------------------------------------*\
                      Class limitedDistanceData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limitedDistanceData
{
    // Private Data

        //- Weight (usually overlap volume)
        List<scalar> weights_;


        // Passive transported data

            //- Max distance
            List<scalar> maxDist_;

            //- Originating point
            List<point> origin_;

            //- Starting data
            List<Type> data_;


    // Private Member Functions


public:


    //- Class used to pass data into container
    class weightFunction
    {
    public:

        const Function1<scalar>& fun_;

        weightFunction(const Function1<scalar>& fun)
        :
            fun_(fun)
        {}

        inline scalar weight
        (
            const polyMesh& mesh,
            const label celli,
            const point& origin,
            const scalar r2
        ) const;
    };


    // Constructors

        //- Construct null with invalid (-1) for distance, null constructor
        //  for data
        inline limitedDistanceData();

        //- Construct from components
        inline limitedDistanceData
        (
            const List<scalar>& weights,
            const List<scalar>& maxDist,
            const List<point>& origin,
            const List<Type>& data
        );


    // Member Functions

        // Access

            inline const List<scalar>& weights() const
            {
                return weights_;
            }
            inline List<scalar>& weights()
            {
                return weights_;
            }

            inline const List<scalar>& maxDist() const
            {
                return maxDist_;
            }
            inline List<scalar>& maxDist()
            {
                return maxDist_;
            }

            inline const List<point>& origin() const
            {
                return origin_;
            }
            inline List<point>& origin()
            {
                return origin_;
            }

            inline const List<Type>& data() const
            {
                return data_;
            }
            inline List<Type>& data()
            {
                return data_;
            }


        // Needed by FaceCellWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const limitedDistanceData<Type>&,
                const scalar,
                TrackingData& td
            ) const;

            //- Convert any absolute coordinates into relative to (patch)face
            //  centre
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Apply rotation matrix to any coordinates
            template<class TrackingData>
            inline void transform
            (
                const polyMesh&,
                const tensor&,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const limitedDistanceData<Type>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const limitedDistanceData<Type>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const limitedDistanceData<Type>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal
            (
                const limitedDistanceData<Type>&,
                TrackingData& td
            ) const;



    // Helper

//        inline scalar calcVolume
//        (
//            const polyMesh& mesh,
//            const label celli,
//            const point& origin,
//            const scalar r2
//        ) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const limitedDistanceData<Type>&) const;

        inline bool operator!=(const limitedDistanceData<Type>&) const;


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const limitedDistanceData<Type>&
        );
        friend Istream& operator>> <Type>
        (
            Istream&,
            limitedDistanceData<Type>&
        );
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

////- Data are contiguous if data type is contiguous
//template<class Type>
//struct is_contiguous<limitedDistanceData<Type>> : is_contiguous<Type> {};
//
////- Data are contiguous label if data type is label
//template<class Type>
//struct is_contiguous_label<limitedDistanceData<Type>> :
//    is_contiguous_label<Type> {};
//
////- Data are contiguous scalar if data type is scalar
//template<class Type>
//struct is_contiguous_scalar<limitedDistanceData<Type>> :
//    is_contiguous_scalar<Type>{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "limitedDistanceData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "limitedDistanceDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
