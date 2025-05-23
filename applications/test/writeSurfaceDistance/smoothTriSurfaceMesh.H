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
    Foam::smoothTriSurfaceMesh

Description
    IOoject and searching on triSurface

    Note: when constructing from dictionary has optional parameters:
        - scale     : scaling factor.
        - tolerance : relative tolerance for doing intersections
                      (see triangle::intersection)
        - minQuality: discard triangles with low quality when getting normal

SourceFiles
    smoothTriSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef smoothTriSurfaceMesh_H
#define smoothTriSurfaceMesh_H

#include "triSurfaceMesh.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class smoothTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class smoothTriSurfaceMesh
:
    public triSurfaceMesh
{
private:

        //- Edges with discontinuity
        PackedBoolList isBorderEdge_;

        //- Points on border edges
        PackedBoolList isPointOnBorderEdge_;

        //- Mark feature edges
        void calcFeatureEdges(const scalar featureAngle);

        //- Calculate local point normal
        vector pointNormal(const label, const label) const;


        //- Disallow default bitwise copy construct
        smoothTriSurfaceMesh(const smoothTriSurfaceMesh&);

        //- Disallow default bitwise assignment
        void operator=(const smoothTriSurfaceMesh&);


public:

    //- Runtime type information
    TypeName("smoothTriSurfaceMesh");


    // Constructors

        //- Construct read, specify feature angle
        smoothTriSurfaceMesh(const IOobject& io, const scalar featureAngle);

        //- Construct from IO and dictionary (used by searchableSurface).
        //  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
        smoothTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~smoothTriSurfaceMesh();


    // Member Functions

        //- From a set of points and indices get the normal
        virtual void getNormal
        (
            const List<pointIndexHit>&,
            vectorField& normal
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
