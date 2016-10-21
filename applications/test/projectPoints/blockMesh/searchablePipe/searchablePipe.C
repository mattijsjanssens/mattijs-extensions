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

#include "searchablePipe.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchablePipe, 0);
    addToRunTimeSelectionTable(searchableSurface, searchablePipe, dict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchablePipe::findNearest
(
    const point& sample,
    const scalar distSqr
) const
{
    pointIndexHit nearCentre = edgeTree_().findNearest(sample, distSqr);

    vector d(sample - nearCentre.hitPoint());

    nearCentre.setPoint(nearCentre.hitPoint() + d/mag(d)*radius_);

    return nearCentre;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchablePipe::searchablePipe
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

Foam::searchablePipe::~searchablePipe()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchablePipe::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


Foam::label Foam::searchablePipe::size() const
{
    return eMeshPtr_().points().size();
}


Foam::tmp<Foam::pointField> Foam::searchablePipe::coordinates() const
{
    return eMeshPtr_().points();
}


void Foam::searchablePipe::boundingSpheres
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


void Foam::searchablePipe::findNearest
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


void Foam::searchablePipe::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchablePipe::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = Zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            // Find nearest on centre line
            pointIndexHit centrePt = edgeTree_().findNearest
            (
                info[i].hitPoint(),
                Foam::magSqr(bounds().span())
            );

            normal[i] = info[i].hitPoint()-centrePt.hitPoint();
            normal[i] /= mag(normal[i]);
        }
    }
}


// ************************************************************************* //
