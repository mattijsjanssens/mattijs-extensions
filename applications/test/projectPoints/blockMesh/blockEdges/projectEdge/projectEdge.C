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

#include "searchableSurfacesQueries.H"
#include "projectEdge.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "pointConstraint.H"
#include "plane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(projectEdge, 0);
    addToRunTimeSelectionTable(blockEdge, projectEdge, Istream);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::projectEdge::findNearest
(
    const pointField& start,
    const scalarField& distSqr,
    pointField& near,
    List<pointConstraint>& constraint
) const
{
    // Multi-surface findNearest

    vectorField normal;
    List<pointIndexHit> info;

    geometry_[surfaces_[0]].findNearest(start, distSqr, info);
    geometry_[surfaces_[0]].getNormal(info, normal);

    // Extract useful info
    near.setSize(info.size());
    forAll(info, i)
    {
        near[i] = info[i].hitPoint();
    }
    constraint.setSize(near.size());

    if (surfaces_.size() == 1)
    {
        constraint = pointConstraint();
        forAll(constraint, i)
        {
            constraint[i].applyConstraint(normal[i]);
        }
    }
    else if (surfaces_.size() >= 2)
    {
        // Work space
        pointField near1;
        vectorField normal1;

        label surfi = 1;
        for (label iter = 0; iter < 10; iter++)
        {
            constraint = pointConstraint();
            forAll(constraint, i)
            {
                constraint[i].applyConstraint(normal[i]);
            }

            // Find intersection with next surface
            const searchableSurface& s = geometry_[surfaces_[surfi]];
            s.findNearest(near, distSqr, info);
            s.getNormal(info, normal1);
            near1.setSize(info.size());
            forAll(info, i)
            {
                near1[i] = info[i].hitPoint();
            }

            // Move to intersection
            forAll(near, pointi)
            {
                plane pl0(near[pointi], normal[pointi]);
                plane pl1(near1[pointi], normal1[pointi]);
                plane::ray r(pl0.planeIntersect(pl1));
                vector n = r.dir() / mag(r.dir());

                vector d(r.refPoint()-near[pointi]);
                d -= (d&n)*n;


                Pout<< "point:" << pointi << " from " << near[pointi];

                near[pointi] += d;
                normal[pointi] = normal1[pointi];

                Pout<< " to:" << near[pointi] << endl;

                //Pout<< "point:" << pointi << endl;
                //Pout<< "    pc was:" << constraint[pointi] << endl;
                //Pout<< "    adding:" << normal1[pointi] << endl;
                constraint[pointi].applyConstraint(normal1[pointi]);
                //Pout<< "    pc now:" << constraint[pointi] << endl;
            }

            // Step to next surface
            surfi = surfaces_.fcIndex(surfi);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::projectEdge::projectEdge
(
    const searchableSurfaces& geometry,
    const pointField& points,
    Istream& is
)
:
    blockEdge(points, is),
    geometry_(geometry)
{
    wordList names(is);
    surfaces_.setSize(names.size());
    forAll(names, i)
    {
        surfaces_[i] = geometry_.findSurfaceID(names[i]);

        if (surfaces_[i] == -1)
        {
            FatalIOErrorInFunction(is)
                << "Cannot find surface " << names[i] << " in geometry"
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::projectEdge::position(const scalar lambda) const
{
    // Initial guess
    const pointField start
    (
        1,
        points_[start_] + lambda * (points_[end_] - points_[start_])
    );

Pout<< "lambda:" << lambda
    << " start:" << start[0] << endl;


    pointField boundaryNear(start);

    if (lambda >= SMALL && lambda < 1.0-SMALL && surfaces_.size())
    {
        List<pointConstraint> boundaryConstraint;
        searchableSurfacesQueries::findNearest
        (
            geometry_,
            surfaces_,
            start,
            scalarField(start.size(), magSqr(points_[end_] - points_[start_])),
            boundaryNear,
            boundaryConstraint
        );
    }

Pout<< "lambda:" << lambda
    << " start:" << start[0]
    << " boundaryNear:" << boundaryNear[0] << endl;

    return boundaryNear[0];
}


// ************************************************************************* //
