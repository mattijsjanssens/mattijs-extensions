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

Class
    Foam::intersectEdge

Description
    Defines the edge from the intersection of two surfaces

SourceFiles
    intersectEdge.C

\*---------------------------------------------------------------------------*/

#ifndef intersectEdge_H
#define intersectEdge_H

#include "blockEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointConstraint;

/*---------------------------------------------------------------------------*\
                           Class intersectEdge Declaration
\*---------------------------------------------------------------------------*/

class intersectEdge
:
    public blockEdge
{
    // Private data

        const searchableSurfaces& geometry_;

        //- The indices of surfaces onto which the points are projected
        labelList surfaces_;
        

//        point p1_, p2_, p3_;
//        cylindricalCS cs_;
//        scalar angle_;
//        scalar radius_;

    // Private Member Functions

        void findNearest
        (
            const pointField& start,
            const scalarField& distSqr,
            pointField& near,
            List<pointConstraint>& constraint
        ) const;

//        //- Calculate the coordinate system, angle and radius
//        cylindricalCS calcAngle();

        //- Disallow default bitwise copy construct
        intersectEdge(const intersectEdge&);

        //- Disallow default bitwise assignment
        void operator=(const intersectEdge&);


public:

    //- Runtime type information
    TypeName("intersect");


    // Constructors

//        //- Construct from components
//        intersectEdge
//        (
//            const pointField& points,
//            const label start, const label end,
//            const point& pMid
//        );

        //- Construct from Istream setting pointsList
        intersectEdge
        (
            const searchableSurfaces& geometry,
            const pointField& points,
            Istream&
        );


    //- Destructor
    virtual ~intersectEdge()
    {}


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const;

        //- Return the length of the curve
        virtual scalar length() const
        {
            NotImplemented;
            return 1;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
