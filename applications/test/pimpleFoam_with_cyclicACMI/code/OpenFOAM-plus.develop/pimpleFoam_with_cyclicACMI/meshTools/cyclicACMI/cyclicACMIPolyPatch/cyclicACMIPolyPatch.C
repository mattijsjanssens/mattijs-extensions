/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "cyclicACMIPolyPatch.H"
#include "SubField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, dictionary);
}

const Foam::scalar Foam::cyclicACMIPolyPatch::tolerance_ = 1e-10;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::cyclicACMIPolyPatch::updateAreas() const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    bool updated = false;

    //Pout<< "cyclicACMIPolyPatch::updateAreas() :"
    //    << " AMITime_:" << AMITime_.eventNo()
    //    << " uptodate:" << mesh.upToDatePoints(AMITime_)
    //    << " mesh.time().timeIndex():" << mesh.time().timeIndex()
    //    << " prevTimeIndex_:" << prevTimeIndex_
    //    << endl;


    // Check if underlying AMI up to date
    if (!mesh.upToDatePoints(AMITime_))
    {
        // This should not happen normally since resetAMI is triggered
        // by any point motion.
        resetAMI();

        updated = true;
    }


    // Check if scaling enabled (and necessary)
    if
    (
        srcScalePtr_.valid()
     && (updated || prevTimeIndex_ != mesh.time().timeIndex())
    )
    {
        const scalar t = boundaryMesh().mesh().time().timeOutputValue();

        // Note: ideally preserve src/tgtMask before clipping to tolerance ...

        srcScaledMask_ =
            min
            (
                scalar(1) - tolerance_,
                max(tolerance_, srcScalePtr_->value(t)*srcMask_)
            );

        if (!tgtScalePtr_.valid())
        {
            tgtScalePtr_= srcScalePtr_.clone(neighbPatch());
        }
        tgtScaledMask_ =
            min
            (
                scalar(1) - tolerance_,
                max(tolerance_, tgtScalePtr_->value(t)*tgtMask_)
            );

        if (debug)
        {
            Pout<< "cyclicACMIPolyPatch::updateAreas : scaling masks"
                << " for " << name() << " mask " << gAverage(srcScaledMask_)
                << " and " << nonOverlapPatch().name()
                << " mask " << gAverage(srcScaledMask_) << endl;
        }

        // Calculate areas from the masks
        updateArea(*this, srcScaledMask_, thisSf_, thisNoSf_);
        updateArea(neighbPatch(), tgtScaledMask_, nbrSf_, nbrNoSf_);

        prevTimeIndex_ = mesh.time().timeIndex();
        AMITime_.setUpToDate();
        updated = true;
    }

    return updated;
}


bool Foam::cyclicACMIPolyPatch::upToDate(const regIOobject& io) const
{
    // Is io up to date with
    // - underlying AMI
    // - scaling
    return io.upToDate(AMITime_);
}


void Foam::cyclicACMIPolyPatch::setUpToDate(regIOobject& io) const
{
    io.setUpToDate();
}


void Foam::cyclicACMIPolyPatch::updateArea
(
    const cyclicACMIPolyPatch& pp,
    const scalarField& mask,
    const vectorField& faceArea,
    const vectorField& noFaceArea
) const
{
    if (mask.size())
    {
        vectorField::subField Sf = pp.faceAreas();
        vectorField::subField noSf = pp.nonOverlapPatch().faceAreas();

        forAll(Sf, facei)
        {
            Sf[facei] = faceArea[facei]*mask[facei];
            noSf[facei] = noFaceArea[facei]*(1.0 - mask[facei]);
        }
    }
}


void Foam::cyclicACMIPolyPatch::normaliseWeights
(
    scalarListList& srcWeights,
    scalarField& srcWeightsSum
) const
{
    forAll(srcWeights, i)
    {
        scalarList& wghts = srcWeights[i];
        if (wghts.size())
        {
            scalar& sum = srcWeightsSum[i];

            forAll(wghts, j)
            {
                wghts[j] /= sum;
            }
            sum = 1.0;
        }
    }
}


void Foam::cyclicACMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod&
) const
{
    if (owner())
    {
        const polyPatch& nonOverlapPatch = this->nonOverlapPatch();
        const cyclicACMIPolyPatch& cp =
            refCast<const cyclicACMIPolyPatch>(this->neighbPatch());

        if (debug)
        {
            Pout<< "cyclicACMIPolyPatch::resetAMI : recalculating weights"
                << " for " << name() << " and " << nonOverlapPatch.name()
                << endl;
        }

        if (boundaryMesh().mesh().hasCellCentres())
        {
            if (debug)
            {
                Pout<< "cyclicACMIPolyPatch::resetAMI : clearing cellCentres"
                    << " for " << name() << " and " << nonOverlapPatch.name()
                    << endl;
            }

            //WarningInFunction
            //    << "The mesh already has cellCentres calculated when"
            //    << " resetting ACMI " << name() << "." << endl
            //    << "This is a problem since ACMI adapts the face areas"
            //    << " (to close cells) so this has" << endl
            //    << "to be done before cell centre calculation." << endl
            //    << "This can happen if e.g. the cyclicACMI is after"
            //    << " any processor patches in the boundary." << endl;
            const_cast<polyMesh&>
            (
                boundaryMesh().mesh()
            ).primitiveMesh::clearGeom();
        }


        // Trigger re-building of faceAreas
        (void)boundaryMesh().mesh().faceAreas();


        // Calculate the AMI using partial face-area-weighted. This leaves
        // the weights as fractions of local areas (sum(weights) = 1 means
        // face is fully covered)
        cyclicAMIPolyPatch::resetAMI
        (
            AMIPatchToPatchInterpolation::imPartialFaceAreaWeight
        );

        AMIPatchToPatchInterpolation& AMI =
            const_cast<AMIPatchToPatchInterpolation&>(this->AMI());

        // Output some stats. AMIInterpolation will have already output the
        // average weights ("sum(weights) min:1 max:1 average:1")
        {
            const scalarField& wghtsSum = AMI.srcWeightsSum();

            label nUncovered = 0;
            label nCovered = 0;
            forAll(wghtsSum, facei)
            {
                scalar sum = wghtsSum[facei];
                if (sum < tolerance_)
                {
                    nUncovered++;
                }
                else if (sum > scalar(1)-tolerance_)
                {
                    nCovered++;
                }
            }
            reduce(nUncovered, sumOp<label>());
            reduce(nCovered, sumOp<label>());
            label nTotal = returnReduce(wghtsSum.size(), sumOp<label>());

            Info<< "ACMI: Patch source uncovered/blended/covered = "
                << nUncovered << ", " << nTotal-nUncovered-nCovered
                << ", " << nCovered << endl;
        }
        {
            const scalarField& wghtsSum = AMI.tgtWeightsSum();

            label nUncovered = 0;
            label nCovered = 0;
            forAll(wghtsSum, facei)
            {
                scalar sum = wghtsSum[facei];
                if (sum < tolerance_)
                {
                    nUncovered++;
                }
                else if (sum > scalar(1)-tolerance_)
                {
                    nCovered++;
                }
            }
            reduce(nUncovered, sumOp<label>());
            reduce(nCovered, sumOp<label>());
            label nTotal = returnReduce(wghtsSum.size(), sumOp<label>());

            Info<< "ACMI: Patch target uncovered/blended/covered = "
                << nUncovered << ", " << nTotal-nUncovered-nCovered
                << ", " << nCovered << endl;
        }

        srcMask_ =
            min(scalar(1) - tolerance_, max(tolerance_, AMI.srcWeightsSum()));

        tgtMask_ =
            min(scalar(1) - tolerance_, max(tolerance_, AMI.tgtWeightsSum()));


        if (srcScalePtr_.valid())
        {
            // Save overlap geometry for later scaling

            thisSf_ = faceAreas();
            thisNoSf_ = nonOverlapPatch.faceAreas();
            nbrSf_ = cp.faceAreas();
            nbrNoSf_ = cp.nonOverlapPatch().faceAreas();
        }

        // Calculate areas from the masks
        updateArea
        (
            *this,
            srcMask_,
            this->faceAreas(),
            nonOverlapPatch.faceAreas()
        );

        updateArea
        (
            cp,
            tgtMask_,
            cp.faceAreas(),
            cp.nonOverlapPatch().faceAreas()
        );

        // Re-normalise the weights since the effect of overlap is already
        // accounted for in the area.
        normaliseWeights(AMI.srcWeights(), AMI.srcWeightsSum());
        normaliseWeights(AMI.tgtWeights(), AMI.tgtWeightsSum());

        // Mark current AMI as up to date with points
        boundaryMesh().mesh().setUpToDatePoints(AMITime_);
    }
}


void Foam::cyclicACMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initGeometry : " << name() << endl;
    }

    // Note: calculates transformation and triggers face centre calculation
    cyclicAMIPolyPatch::initGeometry(pBufs);

    // Initialise the AMI early to make sure we adapt the face areas before the
    // cell centre calculation gets triggered.
    resetAMI();
}


void Foam::cyclicACMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::calcGeometry : " << name() << endl;
    }
    cyclicAMIPolyPatch::calcGeometry(pBufs);
}


void Foam::cyclicACMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initMovePoints : " << name() << endl;
    }

    // Note: calculates transformation and triggers face centre calculation
    cyclicAMIPolyPatch::initMovePoints(pBufs, p);

    // Initialise the AMI early. See initGeometry.
    resetAMI();
}


void Foam::cyclicACMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::movePoints : " << name() << endl;
    }
    cyclicAMIPolyPatch::movePoints(pBufs, p);
}


void Foam::cyclicACMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initUpdateMesh : " << name() << endl;
    }
    cyclicAMIPolyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::updateMesh : " << name() << endl;
    }
    cyclicAMIPolyPatch::updateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::clearGeom()
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::clearGeom : " << name() << endl;
    }
    cyclicAMIPolyPatch::clearGeom();
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::srcMask() const
{
    if (srcScalePtr_.valid())
    {
        return srcScaledMask_;
    }
    else
    {
        return srcMask_;
    }
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::tgtMask() const
{
    if (tgtScalePtr_.valid())
    {
        return tgtScaledMask_;
    }
    else
    {
        return tgtMask_;
    }
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    cyclicAMIPolyPatch(name, size, start, index, bm, patchType, transform),
    nonOverlapPatchName_(word::null),
    nonOverlapPatchID_(-1),
    srcMask_(),
//    (
//        IOobject
//        (
//            "srcMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    tgtMask_(),
//    (
//        IOobject
//        (
//            "tgtMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIRequireMatch_ = false;

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch(name, dict, index, bm, patchType),
    nonOverlapPatchName_(dict.lookup("nonOverlapPatch")),
    nonOverlapPatchID_(-1),
    srcMask_(),
//    (
//        IOobject
//        (
//            "srcMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    tgtMask_(),
//    (
//        IOobject
//        (
//            "tgtMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    srcScalePtr_
    (
        dict.found("scale")
      ? PatchFunction1<scalar>::New(*this, "scale", dict)
      : nullptr
    ),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIRequireMatch_ = false;

    if (nonOverlapPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
//    (
//        IOobject
//        (
//            "srcMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    tgtMask_(),
//    (
//        IOobject
//        (
//            "tgtMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    srcScalePtr_
    (
        pp.srcScalePtr_.valid()
      ? pp.srcScalePtr_.clone(*this)
      : nullptr
    ),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIRequireMatch_ = false;

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& nonOverlapPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    nonOverlapPatchName_(nonOverlapPatchName),
    nonOverlapPatchID_(-1),
    srcMask_(),
//    (
//        IOobject
//        (
//            "srcMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    tgtMask_(),
//    (
//        IOobject
//        (
//            "tgtMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    srcScalePtr_
    (
        pp.srcScalePtr_.valid()
      ? pp.srcScalePtr_.clone(*this)
      : nullptr
    ),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIRequireMatch_ = false;

    if (nonOverlapPatchName_ == name())
    {
        FatalErrorInFunction
            << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
//    (
//        IOobject
//        (
//            "srcMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    tgtMask_(),
//    (
//        IOobject
//        (
//            "tgtMask",
//            boundaryMesh().mesh().pointsInstance(),
//            boundaryMesh().mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        0
//    ),
    srcScalePtr_
    (
        pp.srcScalePtr_.valid()
      ? pp.srcScalePtr_.clone(*this)
      : nullptr
    ),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIRequireMatch_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::~cyclicACMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cyclicACMIPolyPatch& Foam::cyclicACMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];

    // Bit of checking now we know neighbour patch
    if (!owner() && srcScalePtr_.valid())
    {
        WarningInFunction
            << "Ignoring \"scale\" setting in slave patch " << name()
            << endl;
        srcScalePtr_.clear();
        tgtScalePtr_.clear();
    }

    return refCast<const cyclicACMIPolyPatch>(pp);
}


Foam::label Foam::cyclicACMIPolyPatch::nonOverlapPatchID() const
{
    if (nonOverlapPatchID_ == -1)
    {
        nonOverlapPatchID_ =
            this->boundaryMesh().findPatchID(nonOverlapPatchName_);

        if (nonOverlapPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal non-overlapping patch name " << nonOverlapPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        if (nonOverlapPatchID_ < index())
        {
            FatalErrorInFunction
                << "Boundary ordering error: " << type()
                << " patch must be defined prior to its non-overlapping patch"
                << nl
                << type() << " patch: " << name() << ", ID:" << index() << nl
                << "Non-overlap patch: " << nonOverlapPatchName_
                << ", ID:" << nonOverlapPatchID_ << nl
                << exit(FatalError);
        }

        const polyPatch& noPp = this->boundaryMesh()[nonOverlapPatchID_];

        bool ok = true;

        if (size() == noPp.size())
        {
            const scalarField magSf(mag(faceAreas()));
            const scalarField noMagSf(mag(noPp.faceAreas()));

            forAll(magSf, facei)
            {
                scalar ratio = mag(magSf[facei]/(noMagSf[facei] + ROOTVSMALL));

                if (ratio - 1 > tolerance_)
                {
                    ok = false;
                    break;
                }
            }
        }
        else
        {
            ok = false;
        }

        if (!ok)
        {
            FatalErrorInFunction
                << "Inconsistent ACMI patches " << name() << " and "
                << noPp.name() << ".  Patches should have identical topology"
                << exit(FatalError);
        }
    }

    return nonOverlapPatchID_;
}


void Foam::cyclicACMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    cyclicAMIPolyPatch::initOrder(pBufs, pp);
}


bool Foam::cyclicACMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    return cyclicAMIPolyPatch::order(pBufs, pp, faceMap, rotation);
}


void Foam::cyclicACMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    os.writeEntry("nonOverlapPatch", nonOverlapPatchName_);

    if (owner() && srcScalePtr_.valid())
    {
        srcScalePtr_->writeData(os);
    }
}


// ************************************************************************* //
