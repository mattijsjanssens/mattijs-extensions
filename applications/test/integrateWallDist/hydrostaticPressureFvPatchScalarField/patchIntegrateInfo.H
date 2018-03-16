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

Class
    Foam::patchIntegrateInfo

Description

SourceFiles
    patchIntegrateInfoI.H
    patchIntegrateInfo.C

\*---------------------------------------------------------------------------*/

#ifndef patchIntegrateInfo_H
#define patchIntegrateInfo_H

//#include "point.H"
//#include "label.H"
//#include "scalar.H"
//#include "tensor.H"
//#include "pTraits.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class polyMesh;


// Forward declaration of friend functions and operators

class patchIntegrateInfo;

Istream& operator>>(Istream&, patchIntegrateInfo&);
Ostream& operator<<(Ostream&, const patchIntegrateInfo&);


/*---------------------------------------------------------------------------*\
                           Class patchIntegrateInfo Declaration
\*---------------------------------------------------------------------------*/

class patchIntegrateInfo
{
public:

    //- Class used to pass non-cell data to the update function
    class trackData
    {
    public:

        //- Integration distance
        const scalarField& faceCentres_;
        const scalarField& edgeCentres_;

        //- Weight field for integration
        const scalarField& fFld_;
        const scalarField& eFld_;

        // Constructors

            trackData
            (
                const scalarField& faceCentres,
                const scalarField& edgeCentres,
                const scalarField& fFld,
                const scalarField& eFld
            )
            :
                faceCentres_(faceCentres),
                edgeCentres_(edgeCentres),
                fFld_(fFld),
                eFld_(eFld)
            {}
    };


private:

    // Private data

        //- Integrated value
        scalar value_;

public:

    // Constructors

        //- Construct null
        inline patchIntegrateInfo();

        //- Construct from components
        inline patchIntegrateInfo(const scalar);


    // Member Functions

        // Access

            inline scalar value() const;


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Apply rotation matrix
            template<class TrackingData>
            inline void transform
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const label edgeI,
                const label facei,
                const patchIntegrateInfo& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const patchIntegrateInfo& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const primitivePatch& patch,
                const label facei,
                const label edgeI,
                const patchIntegrateInfo& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const patchIntegrateInfo&, TrackingData&) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const patchIntegrateInfo&) const;
        inline bool operator!=(const patchIntegrateInfo&) const;


    // IOstream Operators

        friend inline Ostream& operator<<(Ostream&, const patchIntegrateInfo&);
        friend inline Istream& operator>>(Istream&, patchIntegrateInfo&);
};


//- Data associated with patchIntegrateInfo type are contiguous
template<>
inline bool contiguous<patchIntegrateInfo>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchIntegrateInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
