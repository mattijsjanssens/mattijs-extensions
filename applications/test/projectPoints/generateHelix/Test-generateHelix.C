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

Application

Description

\*---------------------------------------------------------------------------*/

#include "edgeMesh.H"
#include "DynamicField.H"
#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    const point start(-3.5, 0.4, 0.0);

    const scalar radius(0.5);
    const scalar slope(0.2);
    // Number of revolutions (2pi)
    const scalar nRevs(1);

    const scalar tEnd = nRevs*constant::mathematical::twoPi;

    DynamicField<point> points;
    DynamicList<edge> edges;


    // t=0 point
    scalar t = 0.0;
    points.append(point(radius*Foam::sin(t), slope*t, radius*Foam::cos(t)));

    for (; t < tEnd; t += tEnd/100.0)
    {
        label sz = points.size();
        points.append(point(radius*Foam::sin(t), slope*t, radius*Foam::cos(t)));
        edges.append(edge(sz-1, sz));
    }
    edgeMesh eMesh(points+start, edges);

    eMesh.write("spiral.obj");

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
