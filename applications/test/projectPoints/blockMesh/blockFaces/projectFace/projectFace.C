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

#include "projectFace.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "blockDescriptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockFaces
{
    defineTypeNameAndDebug(projectFace, 0);
    addToRunTimeSelectionTable(blockFace, projectFace, Istream);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::searchableSurface& Foam::blockFaces::projectFace::lookupSurface
(
    const searchableSurfaces& geometry,
    Istream& is
) const
{
    word name(is);

    forAll(geometry, i)
    {
        if (geometry[i].name() == name)
        {
            return geometry[i];
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot find surface " << name << " in geometry"
        << exit(FatalIOError);

    return geometry[0];
}


Foam::label Foam::blockFaces::projectFace::index
(
    const labelPair& n,
    const labelPair& coord
) const
{
    return coord.first()+coord.second()*n.first();
}


void Foam::blockFaces::projectFace::calcLambdas
(
    const labelPair& n,
    const pointField& points,
    scalarField& lambdaI,
    scalarField& lambdaJ
) const
{
    lambdaI.setSize(points.size());
    lambdaI = 0.0;
    lambdaJ.setSize(points.size());
    lambdaJ = 0.0;

    for (label i = 1; i < n.first(); i++)
    {
        Pout<< "Row:" << i << endl;

        for (label j = 1; j < n.second(); j++)
        {
            //Pout<< "    Col:" << j << endl;
            label ij = index(n, labelPair(i, j));

            label iMin1j = index(n, labelPair(i-1, j));
            label ijMin1 = index(n, labelPair(i, j-1));

            Pout<< "    ij:" << ij << " iMin1j:" << iMin1j
                << " ijMin1:" << ijMin1
                << " pt:" << points[ij] << endl;

            lambdaI[ij] = lambdaI[iMin1j] + mag(points[ij]-points[iMin1j]);
            lambdaJ[ij] = lambdaJ[ijMin1] + mag(points[ij]-points[ijMin1]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockFaces::projectFace::projectFace
(
    const searchableSurfaces& geometry,
    Istream& is
)
:
    blockFace(is),
    surface_(lookupSurface(geometry, is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockFaces::projectFace::project
(
    const blockDescriptor& desc,
    const label blockFacei,
    pointField& points
) const
{
DebugVar(blockFacei);
DebugVar(desc.curvedFaces()[blockFacei]);
DebugVar(points);
DebugVar(desc.density());
const cellShape& hex = desc.blockShape();
DebugVar(hex);
DebugVar(desc.expand());
DebugVar(desc.vertices());
//DebugVar(desc.facePoints(desc.vertices()));
DebugVar(hex.faces()[blockFacei]);

    label ni = -1;
    label nj = -1;
    switch (blockFacei)
    {
        case 0:
        case 1:
        {
            ni = desc.density()[1]+1;
            nj = desc.density()[2]+1;
        }
        break;

        case 2:
        case 3:
        {
            ni = desc.density()[0]+1;
            nj = desc.density()[2]+1;
        }
        break;

        case 4:
        case 5:
        {
            ni = desc.density()[0]+1;
            nj = desc.density()[1]+1;
        }
        break;
    }

    DebugVar(ni);
    DebugVar(nj);

    const labelPair n(ni, nj);


    scalarField lambdaI(points.size(), 0.0);
    scalarField lambdaJ(points.size(), 0.0);
    calcLambdas(n, points, lambdaI, lambdaJ);
    DebugVar(lambdaI);
    DebugVar(lambdaJ);


    const label nIter = 3;

    for (label iter = 0; iter < nIter;  iter++)
    {
        List<pointIndexHit> hits;
        scalarField nearestDistSqr
        (
            points.size(),
            magSqr(points[0] - points[points.size()-1])
        );
        surface_.findNearest(points, nearestDistSqr, hits);

        forAll(hits, i)
        {
            if (hits[i].hit())
            {
                points[i] = hits[i].hitPoint();
            }
        }

        if (iter < nIter-1)
        {
            scalarField lI;
            scalarField lJ;
            calcLambdas(n, points, lI, lJ);
            DebugVar(lI);
            DebugVar(lJ);

            for (label i = 1; i < n.first()-1; i++)
            {
                for (label j = 1; j < n.second()-1; j++)
                {
                    label ij = index(n, labelPair(i, j));

                    // Predict along i
                    point predi;
                    {
                        label iMin1j = index(n, labelPair(i-1, j));
                        label iLastj = index(n, labelPair(n.first()-1, j));

                        vector v(points[ij]-points[iMin1j]);
                        scalar nearDelta = mag(v)/lI[iLastj];
                        scalar wantedDelta =
                            (lambdaI[ij]-lambdaI[iMin1j])
                           /lambdaI[iLastj];
                        predi = points[iMin1j] + wantedDelta/nearDelta*v;
                    }

                    // Predict along j
                    point predj;
                    {
                        label ijMin1 = index(n, labelPair(i, j-1));
                        label ijLast = index(n, labelPair(i, n.second()-1));

                        vector v(points[ij]-points[ijMin1]);
                        scalar nearDelta = mag(v)/lJ[ijLast];
                        scalar wantedDelta =
                            (lambdaJ[ij]-lambdaJ[ijMin1])
                           /lambdaJ[ijLast];
                        predj = points[ijMin1] + wantedDelta/nearDelta*v;
                    }

                    Pout<< "at i:" << i << " j:" << j
                        << " point was:" << points[ij];
                    points[ij] = 0.5*(predi + predj);
                    Pout<< " point now:" << points[ij] << endl;
                }
            }
        }
    }
}


// ************************************************************************* //
