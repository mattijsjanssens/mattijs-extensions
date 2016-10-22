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

#include "searchableCurve.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "edgeMesh.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "linearInterpolationWeights.H"
#include "quaternion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableCurve, 0);
    addToRunTimeSelectionTable(searchableSurface, searchableCurve, dict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableCurve::findNearest
(
    const point& sample,
    const scalar distSqr
) const
{
    pointIndexHit curvePt = edgeTree_().findNearest(sample, distSqr);

    vector d(sample - curvePt.hitPoint());

    curvePt.setPoint(curvePt.hitPoint() + d/mag(d)*radius_);

    return curvePt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableCurve::searchableCurve
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    eMeshPtr_
    (
        edgeMesh::New
        (
            IOobject
            (
                dict.lookup("file"),                // name
                io.time().constant(),               // instance
                "triSurface",                       // local
                io.time(),                          // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ).objectPath()
        )
    ),
    radius_(readScalar(dict.lookup("radius")))
{
    const edgeMesh& eMesh = eMeshPtr_();

    const pointField& points = eMesh.points();
    const edgeList& edges = eMesh.edges();
    bounds() = boundBox(points, false);

    vector halfSpan(0.5*bounds().span());
    point ctr(bounds().midpoint());

    bounds().min() = ctr - mag(halfSpan)*vector(1, 1, 1);
    bounds().max() = ctr + mag(halfSpan)*vector(1, 1, 1);

    // Calculate bb of all points
    treeBoundBox bb(bounds());

    // Slightly extended bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    edgeTree_.reset
    (
        new indexedOctree<treeDataEdge>
        (
            treeDataEdge
            (
                false,                  // do not cache bb
                edges,
                points,
                identity(edges.size())
            ),
            bb,     // overall search domain
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableCurve::~searchableCurve()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableCurve::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


Foam::label Foam::searchableCurve::size() const
{
    return eMeshPtr_().points().size();
}


Foam::tmp<Foam::pointField> Foam::searchableCurve::coordinates() const
{
    return eMeshPtr_().points();
}


void Foam::searchableCurve::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres = eMeshPtr_().points();
    radiusSqr.setSize(centres.size());
    radiusSqr = Foam::sqr(radius_);
    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}


void Foam::searchableCurve::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableCurve::findInterpolatedNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    const edgeMesh& mesh = eMeshPtr_();
    const indexedOctree<treeDataEdge>& tree = edgeTree_();
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();



    scalarField lambdas(samples.size());
    {
        lambdas[0] = 0.0;
        for (label i = 1; i < samples.size(); i++)
        {
            lambdas[i] = lambdas[i-1] + mag(samples[i]-samples[i-1]);
        }
        lambdas /= lambdas.last();
    }


    pointField curvePoints(samples);
    scalarField curveLambdas(samples.size());
    vectorField axialVecs(samples.size());


    // Upper limit for number of iterations
    const label maxIter = 10;
    // Residual tolerance
    const scalar relTol = 0.1;
    const scalar absTol = 1e-4;

    scalar initialResidual = 0.0;

    for (label iter = 0; iter < maxIter; iter++)
    {
        // Update curve points by projecting onto the curve:
        // - curvePoints
        // - curveLambdas
        // - axialDirs

        forAll(curvePoints, i)
        {
            pointIndexHit curveInfo = tree.findNearest
            (
                curvePoints[i],
                Foam::magSqr(bounds().span())
            );
            curvePoints[i] = curveInfo.hitPoint();

            axialVecs[i] = edges[curveInfo.index()].vec(points);
            axialVecs[i] /= mag(axialVecs[i]);

            if (i == 0)
            {
                curveLambdas[i] = 0.0;
            }
            else
            {
                curveLambdas[i] =
                    curveLambdas[i-1]
                  + mag(curvePoints[i]-curvePoints[i-1]);
            }
        }
        curveLambdas /= curveLambdas.last();

        // Interpolation engine
        linearInterpolationWeights interpolator(curveLambdas);

        // Compare actual distances and move points (along straight line;
        // not along surface)
        vectorField residual(curvePoints.size(), vector::zero);
        labelList indices;
        scalarField weights;
        for (label i = 1; i < curvePoints.size() - 1; i++)
        {
            interpolator.valueWeights(lambdas[i], indices, weights);

            point predicted = vector::zero;
            forAll(indices, indexi)
            {
                predicted += weights[indexi]*curvePoints[indices[indexi]];
            }
            residual[i] = predicted-curvePoints[i];
        }

        scalar scalarResidual = sum(mag(residual));

        if (debug)
        {
            Pout<< "Iter:" << iter << " initialResidual:" << initialResidual
                << " residual:" << scalarResidual << endl;
        }

        if (scalarResidual < absTol*0.5*lambdas.size())
        {
            break;
        }
        else if (iter == 0)
        {
            initialResidual = scalarResidual;
        }
        else if (scalarResidual/initialResidual < relTol)
        {
            break;
        }

        // ? Underrelaxation so as not to overshoot
        curvePoints += 0.5*residual;
    }


    info.setSize(samples.size());
    info = pointIndexHit();

    // Given the current lambdas interpolate radial direction inbetween
    // endpoints
    quaternion qStart;
    vector radialStart;
    {
        radialStart = samples[0]-curvePoints[0];
        radialStart -= (radialStart&axialVecs[0])*axialVecs[0];
        radialStart /= mag(radialStart);
        qStart = quaternion(radialStart, 0.0);

        info[0] = pointIndexHit(true, samples[0], 0);
    }

    quaternion qEnd;
    {
        vector radialEnd(samples.last()-curvePoints.last());
        radialEnd -= (radialEnd&axialVecs.last())*axialVecs.last();
        radialEnd /= mag(radialEnd);
        qEnd = quaternion(radialEnd, 0.0);

        info.last() = pointIndexHit(true, samples.last(), 0);
    }


    for (label i = 1; i < samples.size()-1; i++)
    {
        quaternion q(slerp(qStart, qEnd, curveLambdas[i]));
        vector radialDir(q.transform(radialStart));
        radialDir /= mag(radialDir);

        info[i] = pointIndexHit(true, curvePoints[i]+radius_*radialDir, 0);
    }
}


void Foam::searchableCurve::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableCurve::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    const edgeMesh& mesh = eMeshPtr_();
    const indexedOctree<treeDataEdge>& tree = edgeTree_();
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    normal.setSize(info.size());
    normal = Zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            // Find nearest on curve
            pointIndexHit curvePt = tree.findNearest
            (
                info[i].hitPoint(),
                Foam::magSqr(bounds().span())
            );

            normal[i] = info[i].hitPoint()-curvePt.hitPoint();

            // Subtract axial direction
            vector axialVec = edges[curvePt.index()].vec(points);
            axialVec /= mag(axialVec);
            normal[i] -= (normal[i]&axialVec)*axialVec;

            normal[i] /= mag(normal[i]);
        }
    }
}


// ************************************************************************* //
