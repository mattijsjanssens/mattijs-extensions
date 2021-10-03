/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "cyclicAMIPolyPatch.H"
#include "transformField.H"
#include "SubField.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "faceAreaIntersect.H"
#include "ops.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::vector Foam::cyclicAMIPolyPatch::findFaceNormalMaxRadius
(
    const pointField& faceCentres
) const
{
    // Determine a face furthest away from the axis

    const vectorField n((faceCentres - rotationCentre_) ^ rotationAxis_);

    const scalarField magRadSqr(magSqr(n));

    label facei = findMax(magRadSqr);

    if (debug)
    {
        Info<< "findFaceMaxRadius(const pointField&) : patch: " << name() << nl
            << "    rotFace  = " << facei << nl
            << "    point    = " << faceCentres[facei] << nl
            << "    distance = " << Foam::sqrt(magRadSqr[facei])
            << endl;
    }

    return n[facei];
}


void Foam::cyclicAMIPolyPatch::calcTransforms
(
    const primitivePatch& half0,
    const pointField& half0Ctrs,
    const vectorField& half0Areas,
    const pointField& half1Ctrs,
    const vectorField& half1Areas
)
{
    if (transform() != neighbPatch(0).transform())
    {
        FatalErrorInFunction
            << "Patch " << name()
            << " has transform type " << transformTypeNames[transform()]
            << ", neighbour patch " << neighbPatchName()
            << " has transform type "
            << neighbPatch(0).transformTypeNames[neighbPatch(0).transform()]
            << exit(FatalError);
    }


    // Calculate transformation tensors

    switch (transform())
    {
        case ROTATIONAL:
        {
            tensor revT = Zero;

            if (rotationAngleDefined_)
            {
                const tensor T(rotationAxis_*rotationAxis_);

                const tensor S
                (
                    0, -rotationAxis_.z(), rotationAxis_.y(),
                    rotationAxis_.z(), 0, -rotationAxis_.x(),
                    -rotationAxis_.y(), rotationAxis_.x(), 0
                );

                const tensor revTPos
                (
                    T
                  + cos(rotationAngle_)*(tensor::I - T)
                  + sin(rotationAngle_)*S
                );

                const tensor revTNeg
                (
                    T
                  + cos(-rotationAngle_)*(tensor::I - T)
                  + sin(-rotationAngle_)*S
                );

                // Check - assume correct angle when difference in face areas
                // is the smallest
                const vector transformedAreaPos = gSum(half1Areas & revTPos);
                const vector transformedAreaNeg = gSum(half1Areas & revTNeg);
                const vector area0 = gSum(half0Areas);
                const scalar magArea0 = mag(area0) + ROOTVSMALL;

                // Areas have opposite sign, so sum should be zero when correct
                // rotation applied
                const scalar errorPos = mag(transformedAreaPos + area0);
                const scalar errorNeg = mag(transformedAreaNeg + area0);

                const scalar normErrorPos = errorPos/magArea0;
                const scalar normErrorNeg = errorNeg/magArea0;

                if (errorPos > errorNeg && normErrorNeg < matchTolerance())
                {
                    revT = revTNeg;
                    rotationAngle_ *= -1;
                }
                else
                {
                    revT = revTPos;
                }

                const scalar areaError = min(normErrorPos, normErrorNeg);

                if (areaError > matchTolerance())
                {
                    WarningInFunction
                        << "Patch areas are not consistent within "
                        << 100*matchTolerance()
                        << " % indicating a possible error in the specified "
                        << "angle of rotation" << nl
                        << "    owner patch     : " << name() << nl
                        << "    neighbour patch : " << neighbPatch(0).name()
                        << nl
                        << "    angle           : "
                        << radToDeg(rotationAngle_) << " deg" << nl
                        << "    area error      : " << 100*areaError << " %"
                        << "    match tolerance : " <<  matchTolerance()
                        << endl;
                }

                if (debug)
                {
                    scalar theta = radToDeg(rotationAngle_);

                    Pout<< "cyclicAMIPolyPatch::calcTransforms: patch:"
                        << name()
                        << " Specified rotation:"
                        << " swept angle: " << theta << " [deg]"
                        << " reverse transform: " << revT
                        << endl;
                }
            }
            else
            {
                point n0 = Zero;
                point n1 = Zero;
                if (half0Ctrs.size())
                {
                    n0 = findFaceNormalMaxRadius(half0Ctrs);
                }
                if (half1Ctrs.size())
                {
                    n1 = -findFaceNormalMaxRadius(half1Ctrs);
                }

                reduce(n0, maxMagSqrOp<point>());
                reduce(n1, maxMagSqrOp<point>());

                n0.normalise();
                n1.normalise();

                // Extended tensor from two local coordinate systems calculated
                // using normal and rotation axis
                const tensor E0
                (
                    rotationAxis_,
                    (n0 ^ rotationAxis_),
                    n0
                );
                const tensor E1
                (
                    rotationAxis_,
                    (-n1 ^ rotationAxis_),
                    -n1
                );
                revT = E1.T() & E0;

                if (debug)
                {
                    scalar theta = radToDeg(acos(-(n0 & n1)));

                    Pout<< "cyclicAMIPolyPatch::calcTransforms: patch:"
                        << name()
                        << " Specified rotation:"
                        << " n0:" << n0 << " n1:" << n1
                        << " swept angle: " << theta << " [deg]"
                        << " reverse transform: " << revT
                        << endl;
                }
            }

            const_cast<tensorField&>(forwardT()) = tensorField(1, revT.T());
            const_cast<tensorField&>(reverseT()) = tensorField(1, revT);
            const_cast<vectorField&>(separation()).setSize(0);
            const_cast<boolList&>(collocated()) = boolList(1, false);

            break;
        }
        case TRANSLATIONAL:
        {
            if (debug)
            {
                Pout<< "cyclicAMIPolyPatch::calcTransforms : patch:" << name()
                    << " Specified translation : " << separationVector_
                    << endl;
            }

            const_cast<tensorField&>(forwardT()).clear();
            const_cast<tensorField&>(reverseT()).clear();
            const_cast<vectorField&>(separation()) = vectorField
            (
                1,
                separationVector_
            );
            const_cast<boolList&>(collocated()) = boolList(1, false);

            break;
        }
        default:
        {
            if (debug)
            {
                Pout<< "patch:" << name()
                    << " Assuming cyclic AMI pairs are colocated" << endl;
            }

            const_cast<tensorField&>(forwardT()).clear();
            const_cast<tensorField&>(reverseT()).clear();
            const_cast<vectorField&>(separation()).setSize(0);
            const_cast<boolList&>(collocated()) = boolList(1, true);

            break;
        }
    }

    if (debug)
    {
        Pout<< "patch: " << name() << nl
            << "    forwardT = " << forwardT() << nl
            << "    reverseT = " << reverseT() << nl
            << "    separation = " << separation() << nl
            << "    collocated = " << collocated() << nl << endl;
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::cyclicAMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{
    const labelList& nbrIds = neighbPatchIDs();

    if (debug)
    {
        Pout<< "cyclicAMIPolyPatch::resetAMI : recalculating AMI"
            << " for " << name() << " and neighbours " << nbrIds
            << endl;
    }

    AMIPtrs_.clear();
    AMIPtrs_.setSize(nbrIds.size());
    //AMIIndices_.clear();
    //neighbIndex_.setSize(nbrIds.size());
    //neighbIndex_ = -1;

    forAll(nbrIds, nbri)
    {
        const auto& nbr = neighbPatch(nbri);

        //Pout<< "** cyclicAMIPolyPatch::resetAMI from patch:" << name()
        //    << " patchi:" << this->index()
        //    << " to neighbour:" << nbr.name()
        //    << " patchi:" << nbr.index() << endl;

        if (owner(nbri))
        {
            pointField nbrPoints
            (
                nbr.boundaryMesh().mesh().points(),
                nbr.meshPoints()
            );

            if (debug)
            {
                const Time& t = boundaryMesh().mesh().time();
                OFstream os(t.path()/name() + "_neighbourPatch-org.obj");
                meshTools::writeOBJ(os, nbr.localFaces(), nbrPoints);
            }

            // Transform neighbour patch to local system
            transformPosition(nbrPoints);
            primitivePatch nbrPatch0
            (
                SubList<face>
                (
                    nbr.localFaces(),
                    nbr.size()
                ),
                nbrPoints
            );

            if (debug)
            {
                const Time& t = boundaryMesh().mesh().time();
                OFstream osN(t.path()/name() + "_neighbourPatch-trans.obj");
                meshTools::writeOBJ(osN, nbrPatch0.localFaces(), nbrPoints);

                OFstream osO(t.path()/name() + "_ownerPatch.obj");
                meshTools::writeOBJ(osO, this->localFaces(), localPoints());
            }

            const boundBox srcBb(localPoints(), true);
            const boundBox tgtBb(nbrPoints, true);

            if (!srcBb.overlaps(tgtBb))
            {
                // Currently AMI methods don't allow no overlap so explicitly
                // check for this
                if (debug)
                {
                    Pout<< "cyclicAMIPolyPatch : " << name()
                        << " deleting AMI to " << nbr.name()
                        << " since srcBb:" << srcBb
                        << " tgtBb:" << tgtBb << nl << endl;
                }
            }
            else
            {
                // Construct/apply AMI interpolation to determine addressing and
                // weights. Note that with multiple nbrs the low weight
                // correction is not applied on a per-AMI basis so disabled
                // here (see cyclicAMIPolyPatch::interpolate instead)

                AMIPtrs_.set
                (
                    nbri,
                    new AMIPatchToPatchInterpolation
                    (
                        *this,
                        nbrPatch0,
                        surfPtr(),
                        faceAreaIntersect::tmMesh,
                        AMIRequireMatch_,
                        AMIMethod,
                        (nbrIds.size() == 1 ? AMILowWeightCorrection_ : -1.0),
                        AMIReverse_
                    )
                );

                if
                (
                    gAverage(AMIPtrs_[nbri].tgtWeightsSum()) < SMALL
                 && gAverage(AMIPtrs_[nbri].srcWeightsSum()) < SMALL
                )
                {
                    AMIPtrs_.release(nbri);
                    //neighbIndex_[nbri] = -1;
                    if (debug)
                    {
                        Pout<< "cyclicAMIPolyPatch : " << name()
                            << " deleting AMI to " << nbr.name()
                            << " since zero weights" << nl
                            << endl;
                    }
                }
                else
                {
                    if (debug)
                    {
                        Pout<< "cyclicAMIPolyPatch : " << name()
                            << " constructed AMI to " << nbr.name()
                            << " with " << nl
                            << "    " << "srcAddress:"
                            << AMIPtrs_[nbri].srcAddress().size() << nl
                            << "    " << "tgAddress :"
                            << AMIPtrs_[nbri].tgtAddress().size() << nl << endl;
                    }

                    // Store connection from me to neighbour
                    //const label nbrIndex = nbr.neighbPatchIDs().find(index());
                    //AMIIndices_.append(nbri);
                    //neighbIndex_[nbri] = nbrIndex;
                    // Store connection from neighbour to me
                    //nbr.AMIIndices_.append(nbrIndex);
                    //nbr.neighbIndex_.setSize(nbr.neighbPatchIDs().size(), -1);
                    //nbr.neighbIndex_[nbrIndex] = nbri;
                }
            }
        }
    }


    //- Determine compact addressing
    AMIIndices_.clear();
    neighbIndex_.setSize(nbrIds.size());
    neighbIndex_ = -1;
    {
        forAll(nbrIds, nbri)
        {
            const auto& nbr = neighbPatch(nbri);
            const label nbrIndex = nbr.neighbPatchIDs().find(this->index());

            if (AMI(nbri).valid())
            {
                neighbIndex_[nbri] = nbrIndex;
                AMIIndices_.append(nbri);
            }
            else
            {
                // Check if slave has valid AMI (happens if master has been
                // visited before)
                if (nbr.AMI(nbrIndex).valid())
                {
                    neighbIndex_[nbri] = nbrIndex;
                    AMIIndices_.append(nbri);
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "cyclicAMIPolyPatch::resetAMI : recalculating AMI"
            << " for " << name() << " and neighbours " << nbrIds
            << " now have AMIIndices:" << AMIIndices_
            << endl;
    }
}


const Foam::scalarField& Foam::cyclicAMIPolyPatch::weightsSum() const
{
    const labelList& validAMIs = AMIIndices();

    if (validAMIs.size() == 0)
    {
        //- No valid AMI. Ok as long as not used.
        return weightsSum_;
    }
    else if (validAMIs.size() == 1)
    {
        // Use built-in weights
        const label nbri = validAMIs[0];
        const auto myAMI(AMI(nbri));
        if (myAMI.valid())
        {
            return myAMI().srcWeightsSum();
        }
        else
        {
            const auto nbrAMI(neighbPatch(nbri).AMI(neighbIndex(nbri)));
            return nbrAMI().tgtWeightsSum();
        }
    }
    else
    {
        // Calculate weightsSum if necessary
        if (weightsSum_.size() != size())
        {
            if (debug)
            {
                Pout<< "cyclicAMIPolyPatch::weightsSum : patch:" << name()
                    << " patchi:" << this->index()
                    << " to neighbours:" << neighbPatchNames()
                    << endl;
            }

            weightsSum_.setSize(size());
            weightsSum_ = Zero;

            for (const label nbri : validAMIs)
            {
                const auto myAMI(AMI(nbri));
                const auto& nbr = neighbPatch(nbri);

                if (myAMI.valid())
                {
                    const scalarField& w = myAMI().srcWeightsSum();
                    if (w.size() != this->size())
                    {
                        FatalErrorInFunction << "patch:" << this->name()
                            << " to nbr:" << nbr.name()
                            << exit(FatalError);
                    }
                    weightsSum_ += w;
                }
                else
                {
                    // Neighbour has the AMI
                    const auto nbrAMI(nbr.AMI(neighbIndex(nbri)));
                    const scalarField& w = nbrAMI().tgtWeightsSum();
                    if (w.size() != this->size())
                    {
                        FatalErrorInFunction << "patch:" << nbr.name()
                            << " to *this:" << this->name() << exit(FatalError);
                    }
                    weightsSum_ += w;
                }
            }
        }
        return weightsSum_;
    }
}


void Foam::cyclicAMIPolyPatch::calcTransforms()
{
    const cyclicAMIPolyPatch& half0 = *this;
    vectorField half0Areas(half0.size());
    forAll(half0, facei)
    {
        half0Areas[facei] = half0[facei].areaNormal(half0.points());
    }

    const cyclicAMIPolyPatch& half1 = cyclicAMIPolyPatch::neighbPatch(0);
    vectorField half1Areas(half1.size());
    forAll(half1, facei)
    {
        half1Areas[facei] = half1[facei].areaNormal(half1.points());
    }

    calcTransforms
    (
        half0,
        half0.faceCentres(),
        half0Areas,
        half1.faceCentres(),
        half1Areas
    );

    if (debug)
    {
        Pout<< "calcTransforms() : patch: " << name() << nl
            << "    forwardT = " << forwardT() << nl
            << "    reverseT = " << reverseT() << nl
            << "    separation = " << separation() << nl
            << "    collocated = " << collocated() << nl << endl;
    }
}


void Foam::cyclicAMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    // The AMI is no longer valid. Leave it up to demand-driven calculation
    AMIPtrs_.clear();
    AMIIndices_.clear();
    weightsSum_.clear();
    neighbSize_ = -1;

    polyPatch::initGeometry(pBufs);

    // Early calculation of transforms so e.g. cyclicACMI can use them.
    // Note: also triggers primitiveMesh face centre. Note that cell
    // centres should -not- be calculated
    // since e.g. cyclicACMI override face areas
    calcTransforms();
}


void Foam::cyclicAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    // All geometry done inside initGeometry
}


void Foam::cyclicAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    // The AMI is no longer valid. Leave it up to demand-driven calculation
    AMIPtrs_.clear();
    AMIIndices_.clear();
    weightsSum_.clear();
    neighbSize_ = -1;

    polyPatch::initMovePoints(pBufs, p);

    // See below. Clear out any local geometry
    primitivePatch::movePoints(p);

    // Early calculation of transforms. See above.
    calcTransforms();
}


void Foam::cyclicAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);

    // All transformation tensors already done in initMovePoints
}


void Foam::cyclicAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    // The AMI is no longer valid. Leave it up to demand-driven calculation
    AMIPtrs_.clear();

    polyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::clearGeom()
{
    AMIPtrs_.clear();
    polyPatch::clearGeom();
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
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
    coupledPolyPatch(name, size, start, index, bm, patchType, transform),
    nbrPatchNames_(0),  //word::null),
    nbrPatchIDs_(0),    //-1),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    rotationAngleDefined_(false),
    rotationAngle_(0.0),
    separationVector_(Zero),
    AMIPtrs_(0),
    AMIMethod_(AMIPatchToPatchInterpolation::imFaceAreaWeight),
    AMIReverse_(false),
    AMIRequireMatch_(true),
    AMILowWeightCorrection_(-1.0),
    neighbSize_(-1),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    nbrPatchNames_
    (
        dict.lookupOrDefault<wordList>
        (
            "neighbourPatches",
            wordList(1, dict.lookupOrDefault<word>("neighbourPatch", ""))
        )
    ),
    coupleGroup_(dict),
    nbrPatchIDs_(0),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    rotationAngleDefined_(false),
    rotationAngle_(0.0),
    separationVector_(Zero),
    AMIPtrs_(0),
    AMIMethod_
    (
        AMIPatchToPatchInterpolation::interpolationMethodNames_
        [
            dict.lookupOrDefault
            (
                "method",
                AMIPatchToPatchInterpolation::interpolationMethodNames_
                [
                    AMIPatchToPatchInterpolation::imFaceAreaWeight
                ]
            )
        ]
    ),
    AMIReverse_(dict.lookupOrDefault("flipNormals", false)),
    AMIRequireMatch_(true),
    AMILowWeightCorrection_(dict.lookupOrDefault("lowWeightCorrection", -1.0)),
    neighbSize_(-1),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    if (nbrPatchNames_.empty() && !coupleGroup_.valid())
    {
        FatalIOErrorInFunction(dict)
            << "No \"neighbourPatches\" or \"coupleGroup\" provided."
            << exit(FatalIOError);
    }

    // Check not both neighbourPatches and neighbourPatch provided
    if (dict.found("neighbourPatches") && dict.found("neighbourPatch"))
    {
        FatalIOErrorInFunction(dict)
            << "Cannot both provide 'neighbourPatch' and 'neighbourPatches'"
            << exit(FatalIOError);
    }

    if (nbrPatchNames_.found(name))
    {
        FatalIOErrorInFunction(dict)
            << "Neighbour patch names " << nbrPatchNames_
            << " cannot contain this patch " << name
            << exit(FatalIOError);
    }

    switch (transform())
    {
        case ROTATIONAL:
        {
            dict.readEntry("rotationAxis", rotationAxis_);
            dict.readEntry("rotationCentre", rotationCentre_);
            if (dict.readIfPresent("rotationAngle", rotationAngle_))
            {
                rotationAngleDefined_ = true;
                rotationAngle_ = degToRad(rotationAngle_);

                if (debug)
                {
                    Info<< "rotationAngle: " << rotationAngle_ << " [rad]"
                        <<  endl;
                }
            }

            scalar magRot = mag(rotationAxis_);
            if (magRot < SMALL)
            {
                FatalIOErrorInFunction(dict)
                    << "Illegal rotationAxis " << rotationAxis_ << endl
                    << "Please supply a non-zero vector."
                    << exit(FatalIOError);
            }
            rotationAxis_ /= magRot;

            break;
        }
        case TRANSLATIONAL:
        {
            dict.readEntry("separationVector", separationVector_);
            break;
        }
        default:
        {
            // No additional info required
        }
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    nbrPatchNames_(pp.nbrPatchNames_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchIDs_(), //-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtrs_(0),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    neighbSize_(-1),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const wordList& nbrPatchNames
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    nbrPatchNames_(nbrPatchNames),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchIDs_(), //(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtrs_(0),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    neighbSize_(-1),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{
    if (nbrPatchNames_.found(name()))
    {
        FatalErrorInFunction
            << "Neighbour patch name " << nbrPatchNames_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    nbrPatchNames_(pp.nbrPatchNames_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchIDs_(0),    //(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtrs_(0),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    neighbSize_(-1),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::~cyclicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::cyclicAMIPolyPatch::neighbPatchIDs() const
{
    if (nbrPatchIDs_.empty())
    {
        if (coupleGroup_.valid())
        {
            nbrPatchIDs_ = coupleGroup_.findOtherPatchIDs(*this);

            // Make sure nbrPatchNames is compatible with nbrPatchIDs
            nbrPatchNames_.setSize(nbrPatchIDs_.size());
            forAll(nbrPatchIDs_, i)
            {
                const label patchi = nbrPatchIDs_[i];
                nbrPatchNames_[i] = this->boundaryMesh()[patchi].name();
            }
        }
        else
        {
            const wordList& nbrNames = neighbPatchNames();

            nbrPatchIDs_.setSize(nbrNames.size());
            forAll(nbrNames, i)
            {
                const word& nbrName = nbrNames[i];

                nbrPatchIDs_[i] = this->boundaryMesh().findPatchID(nbrName);

                if (nbrPatchIDs_[i] == -1)
                {
                    FatalErrorInFunction
                        << "Illegal neighbourPatch name " << nbrName
                        << nl << "Valid patch names are "
                        << this->boundaryMesh().names()
                        << exit(FatalError);
                }

                // Check that it is a cyclic AMI patch
                const cyclicAMIPolyPatch& nbrPatch =
                    refCast<const cyclicAMIPolyPatch>
                    (
                        this->boundaryMesh()[nbrPatchIDs_[i]]
                    );

                if (!nbrPatch.neighbPatchNames().found(name()))
                {
                    WarningInFunction
                        << "Patch " << name()
                        << " specifies neighbour patch " << nbrName
                        << nl << " but that in return specifies "
                        << nbrPatch.neighbPatchNames() << endl;
                }
            }
        }
    }
    return nbrPatchIDs_;
}


const Foam::labelList& Foam::cyclicAMIPolyPatch::AMIIndices() const
{
    if (AMIPtrs_.empty())
    {
        resetAMI(AMIMethod_);
    }
    return AMIIndices_;
}


Foam::label Foam::cyclicAMIPolyPatch::neighbIndex(const label index) const
{
    if (AMIPtrs_.empty())
    {
        resetAMI(AMIMethod_);
    }
    return neighbIndex_[index];
}


Foam::label Foam::cyclicAMIPolyPatch::neighbPatchID() const
{
    return neighbPatchIDs()[0];
}


Foam::label Foam::cyclicAMIPolyPatch::neighbSize() const
{
    if (neighbSize_ == -1)
    {
        neighbSize_ = 0;
        for (const label index : AMIIndices())
        {
            neighbSize_ += neighbPatch(index).size();
        }
    }
    return neighbSize_;
}


bool Foam::cyclicAMIPolyPatch::owner(const label index) const
{
    return this->index() < neighbPatchIDs()[index];
}


const Foam::cyclicAMIPolyPatch& Foam::cyclicAMIPolyPatch::neighbPatch
(
    const label index
) const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchIDs()[index]];
    return refCast<const cyclicAMIPolyPatch>(pp);
}


const Foam::autoPtr<Foam::searchableSurface>&
Foam::cyclicAMIPolyPatch::surfPtr() const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if
    (
       !surfPtr_.valid()
     && surfType != "none"
     && (index() < min(neighbPatchIDs()))
    )
    {
        word surfName(surfDict_.lookupOrDefault("name", name()));

        const polyMesh& mesh = boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


Foam::tmpNrc<Foam::AMIPatchToPatchInterpolation>
Foam::cyclicAMIPolyPatch::AMI(const label index) const
{
    if (AMIPtrs_.empty())
    {
        resetAMI(AMIMethod_);
    }

    if (!AMIPtrs_.set(index))
    {
        return nullptr;
    }
    else
    {
        return AMIPtrs_[index];
    }
}


bool Foam::cyclicAMIPolyPatch::applyLowWeightCorrection() const
{
    return AMILowWeightCorrection_ > 0;
}


void Foam::cyclicAMIPolyPatch::transformPosition(pointField& l) const
{
    if (!parallel())
    {
        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(forwardT(), l - rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(forwardT(), l);
        }
    }
    else if (separated())
    {
        // transformPosition gets called on the receiving side,
        // separation gets calculated on the sending side so subtract

        const vectorField& s = separation();
        if (s.size() == 1)
        {
            forAll(l, i)
            {
                l[i] -= s[0];
            }
        }
        else
        {
            l -= s;
        }
    }
}


void Foam::cyclicAMIPolyPatch::transformPosition
(
    point& l,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            forwardT().size() == 1
          ? forwardT()[0]
          : forwardT()[facei]
        );

        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(T, l - rotationCentre_) + rotationCentre_;
        }
        else
        {
            l = Foam::transform(T, l);
        }
    }
    else if (separated())
    {
        const vector& s =
        (
            separation().size() == 1
          ? separation()[0]
          : separation()[facei]
        );

        l -= s;
    }
}


void Foam::cyclicAMIPolyPatch::reverseTransformPosition
(
    point& l,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            reverseT().size() == 1
          ? reverseT()[0]
          : reverseT()[facei]
        );

        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(T, l - rotationCentre_) + rotationCentre_;
        }
        else
        {
            l = Foam::transform(T, l);
        }
    }
    else if (separated())
    {
        const vector& s =
        (
            separation().size() == 1
          ? separation()[0]
          : separation()[facei]
        );

        l += s;
    }
}


void Foam::cyclicAMIPolyPatch::reverseTransformDirection
(
    vector& d,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            reverseT().size() == 1
          ? reverseT()[0]
          : reverseT()[facei]
        );

        d = Foam::transform(T, d);
    }
}


void Foam::cyclicAMIPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{}


void Foam::cyclicAMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::cyclicAMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


Foam::label Foam::cyclicAMIPolyPatch::pointFace
(
    const label facei,
    const vector& n,
    point& p
) const
{
    point prt(p);
    reverseTransformPosition(prt, facei);

    vector nrt(n);
    reverseTransformDirection(nrt, facei);

    label nbrFacei = -1;

    for (const label i : AMIIndices())
    {
        const auto& nbr = neighbPatch(i);
        const auto myAMI(AMI(i));

        if (myAMI.valid())
        {
            nbrFacei = myAMI().tgtPointFace(*this, nbr, nrt, facei, prt);
        }
        else
        {
            const auto nbrAMI(nbr.AMI(neighbIndex(i)));
            nbrFacei = nbrAMI().srcPointFace(nbr, *this, nrt, facei, prt);
        }

        if (nbrFacei != -1)
        {
            p = prt;
            break;
        }
    }

    return nbrFacei;
}


void Foam::cyclicAMIPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    // Write either coupleGroup or neighbourPatches
    if (coupleGroup_.valid())
    {
        coupleGroup_.write(os);
    }
    if (nbrPatchNames_.size() == 1)
    {
        // For backwards compatibility
        os.writeEntry("neighbourPatch", nbrPatchNames_[0]);
    }
    else if (nbrPatchNames_.size())
    {
        os.writeEntry("neighbourPatches", nbrPatchNames_);
    }

    switch (transform())
    {
        case ROTATIONAL:
        {
            os.writeEntry("rotationAxis", rotationAxis_);
            os.writeEntry("rotationCentre", rotationCentre_);

            if (rotationAngleDefined_)
            {
                os.writeEntry("rotationAngle", radToDeg(rotationAngle_));
            }

            break;
        }
        case TRANSLATIONAL:
        {
            os.writeEntry("separationVector", separationVector_);
            break;
        }
        case NOORDERING:
        {
            break;
        }
        default:
        {
            // No additional info to write
        }
    }

    if (AMIMethod_ != AMIPatchToPatchInterpolation::imFaceAreaWeight)
    {
        os.writeEntry
        (
            "method",
            AMIPatchToPatchInterpolation::interpolationMethodNames_
            [
                AMIMethod_
            ]
        );
    }

    if (AMIReverse_)
    {
        os.writeEntry("flipNormals", AMIReverse_);
    }

    if (AMILowWeightCorrection_ > 0)
    {
        os.writeEntry("lowWeightCorrection", AMILowWeightCorrection_);
    }

    if (!surfDict_.empty())
    {
        surfDict_.writeEntry(surfDict_.dictName(), os);
    }
}


// ************************************************************************* //
