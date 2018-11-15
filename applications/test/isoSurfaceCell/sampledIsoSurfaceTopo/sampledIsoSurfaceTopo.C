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

\*---------------------------------------------------------------------------*/

#include "sampledIsoSurfaceTopo.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "isoSurfaceTopo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledIsoSurfaceTopo, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledIsoSurfaceTopo,
        word,
        isoSurfaceTopo
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledIsoSurfaceTopo::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // no update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    prevTimeIndex_ = fvm.time().timeIndex();

    // Clear derived data
    sampledSurface::clearGeom();

    // Optionally read volScalarField
    autoPtr<volScalarField> readFieldPtr_;

    // 1. see if field in database
    // 2. see if field can be read
    const volScalarField* cellFldPtr = nullptr;
    if (fvm.foundObject<volScalarField>(isoField_))
    {
        if (debug)
        {
            InfoInFunction << "Lookup " << isoField_ << endl;
        }

        cellFldPtr = &fvm.lookupObject<volScalarField>(isoField_);
    }
    else
    {
        // Bit of a hack. Read field and store.

        if (debug)
        {
            InfoInFunction
                << "Reading " << isoField_
                << " from time " <<fvm.time().timeName()
                << endl;
        }

        readFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    isoField_,
                    fvm.time().timeName(),
                    fvm,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fvm
            )
        );

        cellFldPtr = readFieldPtr_.operator->();
    }
    const volScalarField& cellFld = *cellFldPtr;

    tmp<pointScalarField> pointFld
    (
        volPointInterpolation::New(fvm).interpolate(cellFld)
    );

    PtrList<isoSurfaceTopo> isos(isoVals_.size());
    forAll(isos, isoi)
    {
        isos.set
        (
            isoi,
            new isoSurfaceTopo
            (
                fvm,
                cellFld.primitiveField(),
                pointFld().primitiveField(),
                isoVals_[isoi],
                (regularise_ ? isoSurfaceTopo::DIAGCELL : isoSurfaceTopo::NONE)
            )
        );
    }


    MeshedSurface<face>& mySurface = const_cast<sampledIsoSurfaceTopo&>(*this);

    if (isos.size() == 1)
    {
        // Straight transfer
        mySurface.transfer(isos[0]);
        meshCells_ = isos[0].meshCells();
    }
    else
    {
        label nFaces = 0;
        label nPoints = 0;
        forAll(isos, isoi)
        {
            nFaces += isos[isoi].size();
            nPoints += isos[isoi].points().size();
        }

        faceList allFaces(nFaces);
        labelList allMeshCells(nFaces);
        pointField allPoints(nPoints);

        nFaces = 0;
        nPoints = 0;
        forAll(isos, isoi)
        {
            isoSurfaceTopo& iso = isos[isoi];

            SubList<face> subAll(allFaces, iso.size(), nFaces);
            subAll = iso;
            // Offset point indices
            if (nPoints > 0)
            {
                forAll(subAll, i)
                {
                    face& f = subAll[i];
                    forAll(f, fp)
                    {
                        f[fp] += nPoints;
                    }
                }
            }
            SubList<label>(allMeshCells, iso.size(), nFaces) = iso.meshCells();
            nFaces += iso.size();

            const pointField& pts = iso.points();
            SubList<point>(allPoints, pts.size(), nPoints) = pts;
            nPoints += pts.size();

            // Clear out asap
            isos.set(isoi, nullptr);
        }

        if (nFaces != allFaces.size() || nPoints != allPoints.size())
        {
            FatalErrorInFunction << "nFaces:" << nFaces
                << " nPoints:" << nPoints << exit(FatalError);
        }


        surfZoneList allZones(1);
        allZones[0] = surfZone
        (
            "allFaces",
            allFaces.size(),    // size
            0,                  // start
            0                   // index
        );

        // Transfer
        mySurface.reset(allPoints.xfer(), allFaces.xfer(), allZones.xfer());
        meshCells_.transfer(allMeshCells);
    }

    // triangulate uses remapFaces()
    // - this is somewhat less efficient since it recopies the faces
    // that we just created, but we probably don't want to do this
    // too often anyhow.
    if (triangulate_)
    {
        mySurface.triangulate();
    }

    if (debug)
    {
        Pout<< "sampledIsoSurfaceTopo::updateGeometry() : constructed iso:"
            << nl
            << "    regularise     : " << regularise_ << nl
            << "    triangulate    : " << triangulate_ << nl
            << "    isoField       : " << isoField_ << nl;
        if (isoVals_.size() == 1)
        {
            Pout<< "    isoValue       : " << isoVals_[0] << nl;
        }
        else
        {
            Pout<< "    isoValues      : " << isoVals_ << nl;
        }
        Pout<< "    points         : " << points().size() << nl
            << "    faces          : " << faces().size() << nl
            << "    cut cells      : " << meshCells_.size() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledIsoSurfaceTopo::sampledIsoSurfaceTopo
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    isoField_(dict.lookup("isoField")),
    isoVals_
    (
        dict.found("isoValues")
      ? scalarField(dict.lookup("isoValues"))
      : scalarField(1, readScalar(dict.lookup("isoValue")))
    ),
    regularise_(dict.lookupOrDefault("regularise", true)),
    triangulate_(dict.lookupOrDefault("triangulate", false)),
    zoneKey_(keyType::null),
    prevTimeIndex_(-1),
    meshCells_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledIsoSurfaceTopo::~sampledIsoSurfaceTopo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledIsoSurfaceTopo::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledIsoSurfaceTopo::expire()
{
    // Clear derived data
    sampledSurface::clearGeom();
    MeshedSurface<face>::clearGeom();

    // already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledIsoSurfaceTopo::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurfaceTopo::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurfaceTopo::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurfaceTopo::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurfaceTopo::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurfaceTopo::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurfaceTopo::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledIsoSurfaceTopo::print(Ostream& os) const
{
    os  << "sampledIsoSurfaceTopo: " << name() << " :"
        << "  field:" << isoField_;
    if (isoVals_.size() == 1)
    {
        os << "  value:" << isoVals_[0];
    }
    else
    {
        os << "  values:" << isoVals_;
    }
    os  << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
