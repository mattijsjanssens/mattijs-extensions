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

    // Precalculate a distribution along straight edge
    const point& startPt = points_[start_];
    const vector d(points_[end_]-startPt);

    const label nPoints = 20;

    pointField start(nPoints);
    lambda_.setSize(nPoints);

    lambda_[0] = 0.0;
    start[0] = startPt;

    scalar dLambda = 1.0/(nPoints-1);
    for (label i = 1; i < lambda_.size(); i++)
    {
        lambda_[i] = lambda_[i-1] + dLambda;
        start[i] = startPt + lambda_[i]*d;
    }

    boundaryNear_.setSize(nPoints);
    boundaryConstraint_.setSize(nPoints);
    searchableSurfacesQueries::findNearest
    (
        geometry_,
        surfaces_,
        start,
        scalarField(start.size(), magSqr(points_[end_] - points_[start_])),
        boundaryNear_,
        boundaryConstraint_
    );

    // Lambda of projected points
    projectedLambda_.setSize(nPoints);
    projectedLambda_[0] = 0.0;
    for (label i = 1; i < boundaryNear_.size(); i++)
    {
        projectedLambda_[i] =
            projectedLambda_[i-1]
          + mag(boundaryNear_[i]-boundaryNear_[i-1]);
    }
    projectedLambda_ /= projectedLambda_.last();

    DebugVar(projectedLambda_);
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
