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

void Foam::cyclicACMIPolyPatch::writeStats
(
    const scalarField& wghtsSum,
    const word& origin,
    const word& patchName
) const
{
    label nUncovered = 0;
    label nCovered = 0;
    for (const scalar sum : wghtsSum)
    {
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
    const label nTotal = returnReduce(wghtsSum.size(), sumOp<label>());

    Info<< "ACMI: Patch " << origin << " uncovered/blended/covered by "
        << patchName << " = "
        << nUncovered << ", " << nTotal-nUncovered-nCovered
        << ", " << nCovered << endl;
}


void Foam::cyclicACMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod&
) const
{
    const polyPatch& nonOverlapPatch = this->nonOverlapPatch();

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

    //const labelList& nbrIds = neighbPatchIDs();
    //srcMasks_.setSize(nbrIds.size());
    //scalarField srcMaskSum(size(), Zero);
    //tgtMasks_.setSize(nbrIds.size());

//XXXXXXXXXX
    // Print owner stats
    for (const label nbri : AMIIndices())
    {
        const auto myAMI(AMI(nbri));

        if (myAMI.valid())
        {
            const polyPatch& nbr = neighbPatch(nbri);

            // Output some stats. AMIInterpolation will have already output the
            // average weights ("sum(weights) min:1 max:1 average:1")
            writeStats(myAMI().srcWeightsSum(), "source", nbr.name());
            writeStats(myAMI().tgtWeightsSum(), "target", nbr.name());
        }
    }


    // Pre-calculate the mask from the weights. Use the summed weights from
    // the cyclicAMIPolyPatch
    maskSum_ = min
    (
        scalar(1) - tolerance_,
        max(tolerance_, weightsSum())
    );
//XXXXXXXXXX


//    for (const label nbri : AMIIndices())
//    {
//        const auto myAMI(AMI(nbri));
//
//        if (myAMI.valid())
//        {
//            const polyPatch& nbr = neighbPatch(nbri);
//            auto& AMI = const_cast<AMIPatchToPatchInterpolation&>(myAMI());
//
//            // Output some stats. AMIInterpolation will have already output the
//            // average weights ("sum(weights) min:1 max:1 average:1")
//
//            writeStats(AMI.srcWeightsSum(), "source", nbr.name());
//            writeStats(AMI.tgtWeightsSum(), "target", nbr.name());
//
//            scalarField& srcMask = srcMasks_[nbri];
//            scalarField& tgtMask = tgtMasks_[nbri];
//
//            srcMask =
//                min
//                (
//                    scalar(1) - tolerance_,
//                    max(tolerance_, AMI.srcWeightsSum())
//                );
//
//            srcMaskSum += srcMask;
//
//            tgtMask =
//                min
//                (
//                    scalar(1) - tolerance_,
//                    max(tolerance_, AMI.tgtWeightsSum())
//                );
//        }
//    }

    // Adapt owner side areas. Note that in uncoupled situations (e.g.
    // decomposePar) srcMask, tgtMask can be zero size.
    if (AMIIndices().size())
    {
        vectorField::subField Sf = faceAreas();
        vectorField::subField noSf = nonOverlapPatch.faceAreas();

        forAll(Sf, facei)
        {
            Sf[facei] *= maskSum_[facei];
            noSf[facei] *= 1.0 - maskSum_[facei];
        }
    }

DebugVar(weightsSum());

    // Re-normalise the weights since the effect of overlap is already
    // accounted for in the area.
    for (const label nbri : AMIIndices())
    {
        const auto myAMI(AMI(nbri));

        if (myAMI.valid())
        {
            auto& AMI = const_cast<AMIPatchToPatchInterpolation&>(myAMI());
            scalarListList& srcWeights = AMI.srcWeights();
            scalarField& srcWeightsSum = const_cast<scalarField&>(weightsSum());

            Pout<< " for " << name() << " normalising srcWeights "
                << flatOutput(srcWeights) << endl;

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
        else
        {
            const auto nbrAMI(neighbPatch(nbri).AMI(neighbIndex(nbri)));
            auto& AMI = const_cast<AMIPatchToPatchInterpolation&>(nbrAMI());
            scalarListList& tgtWeights = AMI.tgtWeights();
            scalarField& tgtWeightsSum = const_cast<scalarField&>(weightsSum());

            Pout<< " for " << name() << " normalising tgtWeights "
                << flatOutput(tgtWeights) << endl;

            forAll(tgtWeights, i)
            {
                scalarList& wghts = tgtWeights[i];
                if (wghts.size())
                {
                    scalar& sum = tgtWeightsSum[i];

                    forAll(wghts, j)
                    {
                        wghts[j] /= sum;
                    }
                    sum = 1.0;
                }
            }
        }
    }

DebugVar(weightsSum());

//    // Adapt slave side areas
//    for (const label nbri : AMIIndices())
//    {
//        const auto myAMI(AMI(nbri));
//
//        if (myAMI.valid())
//        {
//            auto& AMI =
//                const_cast<AMIPatchToPatchInterpolation&>(myAMI());
//
//            const polyPatch& nbr = neighbPatch(nbri);
//            const scalarField& tgtMask = tgtMasks_[nbri];
//            if (tgtMask.size())
//            {
//                const cyclicACMIPolyPatch& cp =
//                    refCast<const cyclicACMIPolyPatch>(nbr);
//                const polyPatch& pp = cp.nonOverlapPatch();
//
//                vectorField::subField Sf = cp.faceAreas();
//                vectorField::subField noSf = pp.faceAreas();
//
//                forAll(Sf, facei)
//                {
//                    Sf[facei] *= tgtMask[facei];
//                    noSf[facei] *= 1.0 - tgtMask[facei];
//                }
//            }
//
//            // Re-normalise the weights since the effect of overlap is already
//            // accounted for in the area.
//            {
//                scalarListList& srcWeights = AMI.srcWeights();
//                scalarField& srcWeightsSum = AMI.srcWeightsSum();
//                forAll(srcWeights, i)
//                {
//                    scalarList& wghts = srcWeights[i];
//                    if (wghts.size())
//                    {
//                        scalar& sum = srcWeightsSum[i];
//
//                        forAll(wghts, j)
//                        {
//                            wghts[j] /= sum;
//                        }
//                        sum = 1.0;
//                    }
//                }
//            }
//            {
//                scalarListList& tgtWeights = AMI.tgtWeights();
//                scalarField& tgtWeightsSum = AMI.tgtWeightsSum();
//                forAll(tgtWeights, i)
//                {
//                    scalarList& wghts = tgtWeights[i];
//                    if (wghts.size())
//                    {
//                        scalar& sum = tgtWeightsSum[i];
//                        forAll(wghts, j)
//                        {
//                            wghts[j] /= sum;
//                        }
//                        sum = 1.0;
//                    }
//                }
//            }
//        }
//    }

    // Set the updated flag
    updated_ = true;

    DebugVar(faceAreas());
    DebugVar(nonOverlapPatch.faceAreas());
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


//const Foam::scalarField&
//Foam::cyclicACMIPolyPatch::srcMask(const label index) const
//{
//    return srcMasks_[index];
//}
//
//
//const Foam::scalarField&
//Foam::cyclicACMIPolyPatch::tgtMask(const label index) const
//{
//    return tgtMasks_[index];
//}


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
    //srcMasks_(0),
    //tgtMasks_(0),
    updated_(false)
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
    //srcMasks_(0),
    //tgtMasks_(0),
    updated_(false)
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
    //srcMasks_(0),
    //tgtMasks_(0),
    updated_(false)
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
    const wordList& nbrPatchNames,
    const word& nonOverlapPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchNames),
    nonOverlapPatchName_(nonOverlapPatchName),
    nonOverlapPatchID_(-1),
    //srcMasks_(0),
    //tgtMasks_(0),
    updated_(false)
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
    //srcMasks_(0),
    //tgtMasks_(0),
    updated_(false)
{
    AMIRequireMatch_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::~cyclicACMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cyclicACMIPolyPatch& Foam::cyclicACMIPolyPatch::neighbPatch
(
    const label index
) const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchIDs()[index]];
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
}


// ************************************************************************* //
