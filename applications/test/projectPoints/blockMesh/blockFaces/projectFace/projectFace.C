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
#include "OBJstream.H"

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
        for (label j = 1; j < n.second(); j++)
        {
            label ij = index(n, labelPair(i, j));
            label iMin1j = index(n, labelPair(i-1, j));
            lambdaI[ij] = lambdaI[iMin1j] + mag(points[ij]-points[iMin1j]);

            label ijMin1 = index(n, labelPair(i, j-1));
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
    // For debugging to tag the output
    static label fIter = 0;

    autoPtr<OBJstream> debugStr;
    if (debug)
    {
        debugStr.reset
        (
            new OBJstream("projectFace_" + Foam::name(fIter++) + ".obj")
        );
        Info<< "Face:" << blockFacei << " on block:" << desc.blockShape()
            << " with verts:" << desc.vertices()
            << " writing lines from start points"
            << " to projected points to " << debugStr().name() << endl;
    }


    // Find out range of vertices in face
    labelPair n(-1, -1);
    switch (blockFacei)
    {
        case 0:
        case 1:
        {
            n.first() = desc.density()[1]+1;
            n.second() = desc.density()[2]+1;
        }
        break;

        case 2:
        case 3:
        {
            n.first() = desc.density()[0]+1;
            n.second() = desc.density()[2]+1;
        }
        break;

        case 4:
        case 5:
        {
            n.first() = desc.density()[0]+1;
            n.second() = desc.density()[1]+1;
        }
        break;
    }


    // Calculate initial edge lengths (= unnormalised u,v coordinates)
    scalarField lambdaI(points.size(), 0.0);
    scalarField lambdaJ(points.size(), 0.0);
    calcLambdas(n, points, lambdaI, lambdaJ);


    // Upper limit for number of iterations
    const label maxIter = 100;
    // Residual tolerance
    const scalar relTol = 0.2;

    scalar initialResidual = 0.0;
    scalar iResidual = 0.0;
    scalar jResidual = 0.0;

    for (label iter = 0; iter < maxIter;  iter++)
    {
        // Do projection
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
                    const point& hitPt = hits[i].hitPoint();
                    if (debugStr.valid())
                    {
                        debugStr().write(linePointRef(points[i], hitPt));
                    }
                    points[i] = hitPt;
                }
            }
        }

        if (debug)
        {
            Pout<< "Iter:" << iter << " initialResidual:" << initialResidual
                << " iResidual+jResidual:" << iResidual+jResidual << endl;
        }


        if (iter > 0 && (iResidual+jResidual)/initialResidual < relTol)
        {
            break;
        }


        // Calculate actual edge lengths from the projected points
        scalarField lI;
        scalarField lJ;
        calcLambdas(n, points, lI, lJ);

        // Predict along i
        vectorField residual(points.size(), vector::zero);

        for (label j = 1; j < n.second()-1; j++)
        {
            label iMaxj = index(n, labelPair(n.first()-1, j));

            for (label i = 1; i < n.first()-1; i++)
            {
                label ij = index(n, labelPair(i, j));

                // How close are we to wanted (normalised) lambda
                scalar nearLambda = lI[ij]/lI[iMaxj];
                scalar wantedLambda = lambdaI[ij]/lambdaI[iMaxj];
                scalar f = wantedLambda/nearLambda;

                // Move in straight line direction. Approximation
                // for local surface normal. Could use constraint
                // instead.
                label iMin1j = index(n, labelPair(i-1, j));
                const point& leftPt(points[iMin1j]);
                vector predicted = leftPt + f*(points[ij]-leftPt);
                residual[ij] = predicted-points[ij];
            }
        }

        if (debugStr.valid())
        {
            forAll(points, i)
            {
                const linePointRef ln(points[i], points[i]+residual[i]);
                debugStr().write(ln);
            }
        }

        iResidual = sum(mag(residual));

        // Update points before doing j. Note: is this needed? Complicates
        // residual checking.
        points += residual;


        // Predict along j
        residual = vector::zero;
        for (label i = 1; i < n.first()-1; i++)
        {
            label ijMax = index(n, labelPair(i, n.second()-1));
            for (label j = 1; j < n.second()-1; j++)
            {
                label ij = index(n, labelPair(i, j));

                scalar nearLambda = lJ[ij]/lJ[ijMax];
                scalar wantedLambda = lambdaJ[ij]/lambdaJ[ijMax];
                scalar f = wantedLambda/nearLambda;

                // Move in straight line direction. Approximation
                // for local surface normal. Could use constraint
                // instead.
                label ijMin1 = index(n, labelPair(i, j-1));
                const point& leftPt(points[ijMin1]);
                vector predicted = leftPt + f*(points[ij]-leftPt);
                residual[ij] = predicted-points[ij];
            }
        }

        if (debugStr.valid())
        {
            forAll(points, i)
            {
                const linePointRef ln(points[i], points[i]+residual[i]);
                debugStr().write(ln);
            }
        }

        jResidual = sum(mag(residual));

        if (iter == 0)
        {
            initialResidual = iResidual + jResidual;
        }

        points += residual;
    }
}


// ************************************************************************* //
