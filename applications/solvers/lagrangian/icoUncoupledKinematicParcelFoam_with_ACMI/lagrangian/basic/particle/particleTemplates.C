/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "IOPosition.H"

#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "cyclicACMIPolyPatch.H"
#include "ACMIWallPolyPatch.H"
#include "wedgePolyPatch.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackData>
void Foam::particle::prepareForParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    // Convert the face index to be local to the processor patch
    facei_ = patchFace(patchi, facei_);
}


template<class TrackData>
void Foam::particle::correctAfterParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    const coupledPolyPatch& ppp =
        refCast<const coupledPolyPatch>(mesh_.boundaryMesh()[patchi]);

    celli_ = ppp.faceCells()[facei_];

    // Have patch transform the position
    ppp.transformPosition(position_, facei_);

    // Transform the properties
    if (!ppp.parallel())
    {
        const tensor& T =
        (
            ppp.forwardT().size() == 1
          ? ppp.forwardT()[0]
          : ppp.forwardT()[facei_]
        );
        transformProperties(T);
    }
    else if (ppp.separated())
    {
        const vector& s =
        (
            (ppp.separation().size() == 1)
          ? ppp.separation()[0]
          : ppp.separation()[facei_]
        );
        transformProperties(-s);
    }

    tetFacei_ = facei_ + ppp.start();

    // Faces either side of a coupled patch have matched base indices,
    // tetPti is specified relative to the base point, already and
    // opposite circulation directions by design, so if the vertices
    // are:
    // source:
    // face    (a b c d e f)
    // fPtI     0 1 2 3 4 5
    //            +
    // destination:
    // face    (a f e d c b)
    // fPtI     0 1 2 3 4 5
    //                  +
    // where a is the base point of the face are matching , and we
    // have fPtI = 1 on the source processor face, i.e. vertex b, then
    // this because of the face circulation direction change, vertex c
    // is the characterising point on the destination processor face,
    // giving the destination fPtI as:
    //     fPtI_d = f.size() - 1 - fPtI_s = 6 - 1 - 1 = 4
    // This relationship can be verified for other points and sizes of
    // face.

    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

    // Reset the face index for the next tracking operation
    if (stepFraction_ > (1.0 - SMALL))
    {
        stepFraction_ = 1.0;
        facei_ = -1;
    }
    else
    {
        facei_ += ppp.start();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::particle::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    IOobject procIO(c.fieldIOobject("origProcId", IOobject::MUST_READ));

    if (procIO.headerOk())
    {
        IOField<label> origProcId(procIO);
        c.checkFieldIOobject(c, origProcId);
        IOField<label> origId(c.fieldIOobject("origId", IOobject::MUST_READ));
        c.checkFieldIOobject(c, origId);

        label i = 0;
        forAllIter(typename CloudType, c, iter)
        {
            particle& p = iter();

            p.origProc_ = origProcId[i];
            p.origId_ = origId[i];
            i++;
        }
    }
}


template<class CloudType>
void Foam::particle::writeFields(const CloudType& c)
{
    // Write the cloud position file
    IOPosition<CloudType> ioP(c);
    ioP.write();

    label np =  c.size();

    IOField<label> origProc
    (
        c.fieldIOobject("origProcId", IOobject::NO_READ),
        np
    );
    IOField<label> origId(c.fieldIOobject("origId", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename CloudType, c, iter)
    {
        origProc[i] = iter().origProc_;
        origId[i] = iter().origId_;
        i++;
    }

    origProc.write();
    origId.write();
}


template<class TrackData>
Foam::label Foam::particle::track(const vector& endPosition, TrackData& td)
{
    facei_ = -1;

    // Tracks to endPosition or stop on boundary
    while (!onBoundary() && stepFraction_ < 1.0 - SMALL)
    {
        stepFraction_ += trackToFace(endPosition, td)*(1.0 - stepFraction_);
    }

    return facei_;
}


template<class TrackData>
Foam::scalar Foam::particle::trackToFace
(
    const vector& endPosition,
    TrackData& td
)
{
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::particleType particleType;

    cloudType& cloud = td.cloud();


    const faceList& pFaces = mesh_.faces();
    const pointField& pPts = mesh_.points();
    const vectorField& pC = mesh_.cellCentres();

    facei_ = -1;

    // Pout<< "Particle " << origId_ << " " << origProc_
    //     << " Tracking from " << position_
    //     << " to " << endPosition
    //     << endl;

    // Pout<< "stepFraction " << stepFraction_ << nl
    //     << "celli " << celli_ << nl
    //     << "tetFacei " << tetFacei_ << nl
    //     << "tetPti " << tetPti_
    //     << endl;

    scalar trackFraction = 0.0;

    // Minimum tetrahedron decomposition of each cell of the mesh into
    // using the cell centre, base point on face, and further two
    // points on the face.  For each face of n points, there are n - 2
    // tets generated.

    // The points for each tet are organised to match those used in the
    // tetrahedron class, supplying them in the order:
    //     Cc, basePt, pA, pB
    // where:
    //   + Cc is the cell centre;
    //   + basePt is the base point on the face;
    //   + pA and pB are the remaining points on the face, such that
    //     the circulation, {basePt, pA, pB} produces a positive
    //     normal by the right-hand rule.  pA and pB are chosen from
    //     tetPti_ do accomplish this depending if the cell owns the
    //     face, tetPti_ is the vertex that characterises the tet, and
    //     is the first vertex on the tet when circulating around the
    //     face. Therefore, the same tetPti represents the same face
    //     triangle for both the owner and neighbour cell.
    //
    // Each tet has its four triangles represented in the same order:
    // 0) tri joining a tet to the tet across the face in next cell.
    //    This is the triangle opposite Cc.
    // 1) tri joining a tet to the tet that is in the same cell, but
    //    belongs to the face that shares the edge of the current face
    //    that doesn't contain basePt.  This is the triangle opposite
    //    basePt.

    // 2) tri joining a tet to the tet that is in the same cell, but
    //    belongs to the face that shares the tet-edge (basePt - pB).
    //    This may be on the same face, or a different one.  This is
    //    the triangle opposite basePt.  This is the triangle opposite
    //    pA.

    // 4) tri joining a tet to the tet that is in the same cell, but
    //    belongs to the face that shares the tet-edge (basePt - pA).
    //    This may be on the same face, or a different one.  This is
    //    the triangle opposite basePt.  This is the triangle opposite
    //    pA.

    // Which tri (0..3) of the tet has been crossed
    label triI = -1;

    // Determine which face was actually crossed.  lambdaMin < SMALL
    // is considered a trigger for a tracking correction towards the
    // current tet centre.
    scalar lambdaMin = VGREAT;

    DynamicList<label>& tris = cloud.labels();

    // Tet indices that will be set by hitWallFaces if a wall face is
    // to be hit, or are set when any wall tri of a tet is hit.
    // Carries the description of the tet on which the cell face has
    // been hit.  For the case of being set in hitWallFaces, this may
    // be a different tet to the one that the particle occupies.
    tetIndices faceHitTetIs;

    // What tolerance is appropriate the minimum lambda numerator and
    // denominator for tracking in this cell.
    scalar lambdaDistanceTolerance =
        lambdaDistanceToleranceCoeff*mesh_.cellVolumes()[celli_];

    do
    {
        if (triI != -1)
        {
            // Change tet ownership because a tri face has been crossed
            tetNeighbour(triI);
        }

        const Foam::face& f = pFaces[tetFacei_];

        bool own = (mesh_.faceOwner()[tetFacei_] == celli_);

        label tetBasePtI = mesh_.tetBasePtIs()[tetFacei_];

        label basePtI = f[tetBasePtI];

        label facePtI = (tetPti_ + tetBasePtI) % f.size();
        label otherFacePtI = f.fcIndex(facePtI);

        label fPtAI = -1;
        label fPtBI = -1;

        if (own)
        {
            fPtAI = facePtI;
            fPtBI = otherFacePtI;
        }
        else
        {
            fPtAI = otherFacePtI;
            fPtBI = facePtI;
        }

        tetPointRef tet
        (
            pC[celli_],
            pPts[basePtI],
            pPts[f[fPtAI]],
            pPts[f[fPtBI]]
        );

        if (!tet.inside(position_))
        {
            meshTools::writeOBJ(Pout, tet.a());
            meshTools::writeOBJ(Pout, tet.b());
            meshTools::writeOBJ(Pout, tet.c());
            meshTools::writeOBJ(Pout, tet.d());

            Pout<< "f 1 3 2" << nl
                << "f 2 3 4" << nl
                << "f 1 4 3" << nl
                << "f 1 2 4" << endl;
            meshTools::writeOBJ(Pout, position_);
            meshTools::writeOBJ(Pout, endPosition);
        }


        if (lambdaMin < SMALL)
        {
            // Apply tracking correction towards tet centre

            if (debug)
            {
                Pout<< "tracking rescue using tetCentre from " << position();
            }

            position_ += trackingCorrectionTol*(tet.centre() - position_);

            if (debug)
            {
                Pout<< " to " << position() << " due to "
                    << (tet.centre() - position_) << endl;
            }

            cloud.trackingRescue();

            return trackFraction;
        }

        if (triI != -1 && mesh_.moving())
        {
            // Mesh motion requires stepFraction to be correct for
            // each tracking portion, so trackToFace must return after
            // every lambda calculation.
            return trackFraction;
        }

        FixedList<vector, 4> tetAreas;

        tetAreas[0] = tet.Sa();
        tetAreas[1] = tet.Sb();
        tetAreas[2] = tet.Sc();
        tetAreas[3] = tet.Sd();

        FixedList<label, 4> tetPlaneBasePtIs;

        tetPlaneBasePtIs[0] = basePtI;
        tetPlaneBasePtIs[1] = f[fPtAI];
        tetPlaneBasePtIs[2] = basePtI;
        tetPlaneBasePtIs[3] = basePtI;

        findTris
        (
            endPosition,
            tris,
            tet,
            tetAreas,
            tetPlaneBasePtIs,
            lambdaDistanceTolerance
        );

        // Reset variables for new track
        triI = -1;
        lambdaMin = VGREAT;

        // Pout<< "tris " << tris << endl;

        // Sets a value for lambdaMin and facei_ if a wall face is hit
        // by the track.
        hitWallFaces
        (
            cloud,
            position_,
            endPosition,
            lambdaMin,
            faceHitTetIs
        );

        // Did not hit any tet tri faces, and no wall face has been
        // found to hit.
        if (tris.empty() && facei_ < 0)
        {
            position_ = endPosition;

Pout<< "updated position to:" << position_
    << " since reached end" << endl;


            return 1.0;
        }
        else
        {
            // Loop over all found tris and see if any of them find a
            // lambda value smaller than that found for a wall face.
            forAll(tris, i)
            {
                label tI = tris[i];

                scalar lam = tetLambda
                (
                    position_,
                    endPosition,
                    tI,
                    tetAreas[tI],
                    tetPlaneBasePtIs[tI],
                    celli_,
                    tetFacei_,
                    tetPti_,
                    lambdaDistanceTolerance
                );

                if (lam < lambdaMin)
                {
                    lambdaMin = lam;

                    triI = tI;
                }
            }
        }

        if (triI == 0)
        {
            // This must be a cell face crossing
            facei_ = tetFacei_;

            // Set the faceHitTetIs to those for the current tet in case a
            // wall interaction is required with the cell face
            faceHitTetIs = tetIndices
            (
                celli_,
                tetFacei_,
                tetBasePtI,
                fPtAI,
                fPtBI,
                tetPti_
            );
        }
        else if (triI > 0)
        {
            // A tri was found to be crossed before a wall face was hit (if any)
            facei_ = -1;
        }

        Pout<< "track loop " << position_ << " " << endPosition << nl
             << "    " << celli_
             << "    " << facei_
             << " " << tetFacei_
             << " " << tetPti_
             << " " << triI
             << " " << lambdaMin
             << " " << trackFraction
             << endl;

        // Pout<< "# Tracking loop tet "
        //     << origId_ << " " << origProc_<< nl
        //     << "# face: " << tetFacei_ << nl
        //     << "# tetPti: " << tetPti_ << nl
        //     << "# tetBasePtI: " << mesh_.tetBasePtIs()[tetFacei_] << nl
        //     << "# tet.mag(): " << tet.mag() << nl
        //     << "# tet.quality(): " << tet.quality()
        //     << endl;

        // meshTools::writeOBJ(Pout, tet.a());
        // meshTools::writeOBJ(Pout, tet.b());
        // meshTools::writeOBJ(Pout, tet.c());
        // meshTools::writeOBJ(Pout, tet.d());

        // Pout<< "f 1 3 2" << nl
        //     << "f 2 3 4" << nl
        //     << "f 1 4 3" << nl
        //     << "f 1 2 4" << endl;

        // The particle can be 'outside' the tet.  This will yield a
        // lambda larger than 1, or smaller than 0.  For values < 0,
        // the particle travels away from the tet and we don't move
        // the particle, only change tet/cell.  For values larger than
        // 1, we move the particle to endPosition before the tet/cell
        // change.
        if (lambdaMin > SMALL)
        {
            if (lambdaMin <= 1.0)
            {
                trackFraction += lambdaMin*(1 - trackFraction);
                position_ += lambdaMin*(endPosition - position_);

Pout<< "updated position to:" << position_ << " trackFraction:" << trackFraction
    << " since lamddaMin:" << lambdaMin << endl;
            }
            else
            {
                position_ = endPosition;

Pout<< "updated position to:" << position_ << " trackFraction:" << trackFraction
    << " since reached end" << endl;

                return 1.0;
            }
        }
        else
        {
            // Set lambdaMin to zero to force a towards-tet-centre
            // correction.
            lambdaMin = 0.0;
Pout<< "Froce tet-correction" << endl;
        }

    } while (facei_ < 0);

    particleType& p = static_cast<particleType&>(*this);
    p.hitFace(td);

    if (internalFace(facei_))
    {
        // Change tet ownership because a tri face has been crossed,
        // in general this is:
        //     tetNeighbour(triI);
        // but triI must be 0;
        // No modifications are required for triI = 0, no call required to
        //     tetNeighbour(0);

        if (celli_ == mesh_.faceOwner()[facei_])
        {
            celli_ = mesh_.faceNeighbour()[facei_];
        }
        else if (celli_ == mesh_.faceNeighbour()[facei_])
        {
            celli_ = mesh_.faceOwner()[facei_];
        }
        else
        {
            FatalErrorInFunction
                << "addressing failure" << abort(FatalError);
        }
    }
    else
    {
        label origFacei = facei_;
        label patchi = patch(facei_);

        // No action taken for tetPti_ for tetFacei_ here, handled by
        // patch interaction call or later during processor transfer.

        if
        (
            !p.hitPatch
            (
                mesh_.boundaryMesh()[patchi],
                td,
                patchi,
                trackFraction,
                faceHitTetIs
            )
        )
        {
            // Did patch interaction model switch patches?
            if (facei_ != origFacei)
            {
                patchi = patch(facei_);
            }

            const polyPatch& patch = mesh_.boundaryMesh()[patchi];

            if (isA<wedgePolyPatch>(patch))
            {
                p.hitWedgePatch
                (
                    static_cast<const wedgePolyPatch&>(patch), td
                );
            }
            else if (isA<symmetryPlanePolyPatch>(patch))
            {
                p.hitSymmetryPlanePatch
                (
                    static_cast<const symmetryPlanePolyPatch&>(patch), td
                );
            }
            else if (isA<symmetryPolyPatch>(patch))
            {
                p.hitSymmetryPatch
                (
                    static_cast<const symmetryPolyPatch&>(patch), td
                );
            }
            else if (isA<cyclicPolyPatch>(patch))
            {
                p.hitCyclicPatch
                (
                    static_cast<const cyclicPolyPatch&>(patch), td
                );
            }
            else if (isA<cyclicAMIPolyPatch>(patch))
            {
                p.hitCyclicAMIPatch
                (
                    static_cast<const cyclicAMIPolyPatch&>(patch),
                    td,
                    endPosition - position_
                );
            }
            else if (isA<processorPolyPatch>(patch))
            {
                p.hitProcessorPatch
                (
                    static_cast<const processorPolyPatch&>(patch), td
                );
            }
            else if (isA<ACMIWallPolyPatch>(patch))
            {
                p.hitACMIWallPatch
                (
                    static_cast<const ACMIWallPolyPatch&>(patch),
                    td,
                    faceHitTetIs,
                    endPosition - position_
                );
            }
            else if (isA<wallPolyPatch>(patch))
            {
                p.hitWallPatch
                (
                    static_cast<const wallPolyPatch&>(patch), td, faceHitTetIs
                );
            }
            else
            {
                p.hitPatch(patch, td);
            }
        }
    }

    if (lambdaMin < SMALL)
    {
        // Apply tracking correction towards tet centre.
        // Generate current tet to find centre to apply correction.

        tetPointRef tet = currentTet();

        if (debug)
        {
            Pout<< "tracking rescue for lambdaMin:" << lambdaMin
                << "from " << position();
        }

        position_ += trackingCorrectionTol*(tet.centre() - position_);

        if
        (
            cloud.hasWallImpactDistance()
         && !internalFace(faceHitTetIs.face())
         && cloud.cellHasWallFaces()[faceHitTetIs.cell()]
        )
        {
            const polyBoundaryMesh& patches = mesh_.boundaryMesh();

            label fI = faceHitTetIs.face();

            label patchi = patches.patchID()[fI - mesh_.nInternalFaces()];

            if (isA<wallPolyPatch>(patches[patchi]))
            {
                // In the case of collision with a wall where there is
                // a non-zero wallImpactDistance, it is possible for
                // there to be a tracking correction required to bring
                // the particle into the domain, but the position of
                // the particle is further from the wall than the tet
                // centre, in which case the normal correction can be
                // counter-productive, i.e. pushes the particle
                // further out of the domain.  In this case it is the
                // position that hit the wall that is in need of a
                // rescue correction.

                triPointRef wallTri = faceHitTetIs.faceTri(mesh_);

                tetPointRef wallTet = faceHitTetIs.tet(mesh_);

                vector nHat = wallTri.normal();
                nHat /= mag(nHat);

                const scalar r = p.wallImpactDistance(nHat);

                // Removing (approximately) the wallTri normal
                // component of the existing correction, to avoid the
                // situation where the existing correction in the wall
                // normal direction is larger towards the wall than
                // the new correction is away from it.
                position_ +=
                    trackingCorrectionTol
                   *(
                        (wallTet.centre() - (position_ + r*nHat))
                      - (nHat & (tet.centre() - position_))*nHat
                    );
            }
        }

        if (debug)
        {
            Pout<< " to " << position() << endl;
        }

        cloud.trackingRescue();
    }

    return trackFraction;
}


template<class CloudType>
void Foam::particle::hitWallFaces
(
    const CloudType& cloud,
    const vector& from,
    const vector& to,
    scalar& lambdaMin,
    tetIndices& closestTetIs
)
{
    typedef typename CloudType::particleType particleType;

    if (!(cloud.hasWallImpactDistance() && cloud.cellHasWallFaces()[celli_]))
    {
        return;
    }

    particleType& p = static_cast<particleType&>(*this);

    const faceList& pFaces = mesh_.faces();

    const Foam::cell& thisCell = mesh_.cells()[celli_];

    scalar lambdaDistanceTolerance =
        lambdaDistanceToleranceCoeff*mesh_.cellVolumes()[celli_];

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(thisCell, cFI)
    {
        label fI = thisCell[cFI];

        if (internalFace(fI))
        {
            continue;
        }

        label patchi = patches.patchID()[fI - mesh_.nInternalFaces()];

        if (isA<wallPolyPatch>(patches[patchi]))
        {
            // Get the decomposition of this wall face

            const List<tetIndices> faceTetIs =
                polyMeshTetDecomposition::faceTetIndices(mesh_, fI, celli_);

            const Foam::face& f = pFaces[fI];

            forAll(faceTetIs, tI)
            {
                const tetIndices& tetIs = faceTetIs[tI];

                triPointRef tri = tetIs.faceTri(mesh_);

                vector n = tri.normal();

                vector nHat = n/mag(n);

                // Radius of particle with respect to this wall face
                // triangle.  Assuming that the wallImpactDistance
                // does not change as the particle or the mesh moves
                // forward in time.
                scalar r = p.wallImpactDistance(nHat);

                vector toPlusRNHat = to + r*nHat;

                // triI = 0 because it is the cell face tri of the tet
                // we are concerned with.
                scalar tetClambda = tetLambda
                (
                    tetIs.tet(mesh_).centre(),
                    toPlusRNHat,
                    0,
                    n,
                    f[tetIs.faceBasePt()],
                    celli_,
                    fI,
                    tetIs.tetPt(),
                    lambdaDistanceTolerance
                );

                if ((tetClambda <= 0.0) || (tetClambda >= 1.0))
                {
                    // toPlusRNHat is not on the outside of the plane of
                    // the wall face tri, the tri cannot be hit.
                    continue;
                }

                // Check if the actual trajectory of the near-tri
                // points intersects the triangle.

                vector fromPlusRNHat = from + r*nHat;

                // triI = 0 because it is the cell face tri of the tet
                // we are concerned with.
                scalar lambda = tetLambda
                (
                    fromPlusRNHat,
                    toPlusRNHat,
                    0,
                    n,
                    f[tetIs.faceBasePt()],
                    celli_,
                    fI,
                    tetIs.tetPt(),
                    lambdaDistanceTolerance
                );

                pointHit hitInfo(Zero);

                if (mesh_.moving())
                {
                    // For a moving mesh, the position of wall
                    // triangle needs to be moved in time to be
                    // consistent with the moment defined by the
                    // current value of stepFraction and the value of
                    // lambda just calculated.

                    // Total fraction thought the timestep of the
                    // motion, including stepFraction before the
                    // current tracking step and the current
                    // lambda
                    // i.e.
                    // let s = stepFraction, l = lambda
                    // Motion of x in time:
                    // |-----------------|---------|---------|
                    // x00               x0        xi        x
                    //
                    // where xi is the correct value of x at the required
                    // tracking instant.
                    //
                    // x0 = x00 + s*(x - x00) = s*x + (1 - s)*x00
                    //
                    // i.e. the motion covered by previous tracking portions
                    // within this timestep, and
                    //
                    // xi = x0 + l*(x - x0)
                    //    = l*x + (1 - l)*x0
                    //    = l*x + (1 - l)*(s*x + (1 - s)*x00)
                    //    = (s + l - s*l)*x + (1 - (s + l - s*l))*x00
                    //
                    // let m = (s + l - s*l)
                    //
                    // xi = m*x + (1 - m)*x00 = x00 + m*(x - x00);
                    //
                    // In the same form as before.

                    // Clip lambda to 0.0-1.0 to ensure that sensible
                    // positions are used for triangle intersections.
                    scalar lam = max(0.0, min(1.0, lambda));

                    scalar m = stepFraction_ + lam - (stepFraction_*lam);

                    triPointRef tri00 = tetIs.oldFaceTri(mesh_);

                    // Use SMALL positive tolerance to make the triangle
                    // slightly "fat" to improve robustness.  Intersection
                    // is calculated as the ray (from + r*nHat) -> (to +
                    // r*nHat).

                    point tPtA = tri00.a() + m*(tri.a() - tri00.a());
                    point tPtB = tri00.b() + m*(tri.b() - tri00.b());
                    point tPtC = tri00.c() + m*(tri.c() - tri00.c());

                    triPointRef t(tPtA, tPtB, tPtC);

                    // The point fromPlusRNHat + m*(to - from) is on the
                    // plane of the triangle.  Determine the
                    // intersection with this triangle by testing if
                    // this point is inside or outside of the triangle.
                    hitInfo = t.intersection
                    (
                        fromPlusRNHat + m*(to - from),
                        t.normal(),
                        intersection::FULL_RAY,
                        SMALL
                    );
                }
                else
                {
                    // Use SMALL positive tolerance to make the triangle
                    // slightly "fat" to improve robustness.  Intersection
                    // is calculated as the ray (from + r*nHat) -> (to +
                    // r*nHat).
                    hitInfo = tri.intersection
                    (
                        fromPlusRNHat,
                        (to - from),
                        intersection::FULL_RAY,
                        SMALL
                    );
                }

                if (hitInfo.hit())
                {
                    if (lambda < lambdaMin)
                    {
                        lambdaMin = lambda;

                        facei_ = fI;

                        closestTetIs = tetIs;
                    }
                }
            }
        }
    }
}


template<class TrackData>
void Foam::particle::hitFace(TrackData&)
{}


template<class TrackData>
bool Foam::particle::hitPatch
(
    const polyPatch&,
    TrackData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


template<class TrackData>
void Foam::particle::hitWedgePatch
(
    const wedgePolyPatch& wpp,
    TrackData&
)
{
    FatalErrorInFunction
        << "Hitting a wedge patch should not be possible."
        << abort(FatalError);

    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackData>
void Foam::particle::hitSymmetryPlanePatch
(
    const symmetryPlanePolyPatch& spp,
    TrackData&
)
{
    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackData>
void Foam::particle::hitSymmetryPatch
(
    const symmetryPolyPatch& spp,
    TrackData&
)
{
    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackData>
void Foam::particle::hitCyclicPatch
(
    const cyclicPolyPatch& cpp,
    TrackData& td
)
{
    facei_ = cpp.transformGlobalFace(facei_);

    celli_ = mesh_.faceOwner()[facei_];

    tetFacei_ = facei_;

    // See note in correctAfterParallelTransfer for tetPti_ addressing.
    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

    const cyclicPolyPatch& receiveCpp = cpp.neighbPatch();
    label patchFacei = receiveCpp.whichFace(facei_);

    // Now the particle is on the receiving side

    // Have patch transform the position
    receiveCpp.transformPosition(position_, patchFacei);

    // Transform the properties
    if (!receiveCpp.parallel())
    {
        const tensor& T =
        (
            receiveCpp.forwardT().size() == 1
          ? receiveCpp.forwardT()[0]
          : receiveCpp.forwardT()[patchFacei]
        );
        transformProperties(T);
    }
    else if (receiveCpp.separated())
    {
        const vector& s =
        (
            (receiveCpp.separation().size() == 1)
          ? receiveCpp.separation()[0]
          : receiveCpp.separation()[patchFacei]
        );
        transformProperties(-s);
    }
}


template<class TrackData>
void Foam::particle::hitCyclicAMIPatch
(
    const cyclicAMIPolyPatch& cpp,
    TrackData& td,
    const vector& direction
)
{
    DebugVar(cpp.name());

    const cyclicAMIPolyPatch& receiveCpp = cpp.neighbPatch();

    // Patch face index on sending side
    label patchFacei = facei_ - cpp.start();

    // Patch face index on receiving side - also updates position
    patchFacei = cpp.pointFace(patchFacei, direction, position_);

    if (patchFacei < 0)
    {
        // If the patch face of the particle is not known assume that
        // the particle is lost and to be deleted
        td.keepParticle = false;

        WarningInFunction
            << "Particle lost across " << cyclicAMIPolyPatch::typeName
            << " patches " << cpp.name() << " and " << receiveCpp.name()
            << " at position " << position_
            << endl;
    }

    // Convert face index into global numbering
    facei_ = patchFacei + receiveCpp.start();

    celli_ = mesh_.faceOwner()[facei_];

    tetFacei_ = facei_;

    // See note in correctAfterParallelTransfer for tetPti_ addressing.
    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

    // Now the particle is on the receiving side

    // Transform the properties
    if (!receiveCpp.parallel())
    {
        const tensor& T =
        (
            receiveCpp.forwardT().size() == 1
          ? receiveCpp.forwardT()[0]
          : receiveCpp.forwardT()[patchFacei]
        );
        transformProperties(T);
    }
    else if (receiveCpp.separated())
    {
        const vector& s =
        (
            (receiveCpp.separation().size() == 1)
          ? receiveCpp.separation()[0]
          : receiveCpp.separation()[patchFacei]
        );
        transformProperties(-s);
    }
}


template<class TrackData>
void Foam::particle::hitProcessorPatch(const processorPolyPatch&, TrackData&)
{}


//XXXXXX
template<class TrackData>
void Foam::particle::hitACMIWallPatch
(
    const ACMIWallPolyPatch& pp,
    TrackData& td,
    const tetIndices& ti,
    const vector& direction
)
{
    DebugVar(position_);
    DebugVar(pp.name());
    DebugVar(pp.start());
    DebugVar(pp.size());

    const cyclicACMIPolyPatch& cpp = pp.ACMI();
    const cyclicAMIPolyPatch& receiveCpp = cpp.neighbPatch();

    DebugVar(cpp.name());
    DebugVar(cpp.start());
    DebugVar(cpp.size());

    DebugVar(receiveCpp.name());
    DebugVar(receiveCpp.start());
    DebugVar(receiveCpp.size());


    // Patch face index on sending side
    label patchFacei = facei_ - pp.start();

    DebugVar(patchFacei);

    // Patch face index on receiving side - also updates position
    patchFacei = cpp.pointFace(patchFacei, direction, position_);

    DebugVar(patchFacei);

    if (patchFacei < 0)
    {
        Pout<< "** missed **" << endl;
        DebugVar(*this);
        hitWallPatch(pp, td, ti);
    }
    else
    {
        // Convert face index into global numbering
        facei_ = patchFacei + receiveCpp.start();

        celli_ = mesh_.faceOwner()[facei_];

DebugVar(celli_);
Pout<< "Moved particle through " << receiveCpp.name()
    << " into cell:" << celli_
    << " at:" << mesh_.cellCentres()[celli_]
    << endl;
Pout<< "particle now at:" << position_ << endl;

        tetFacei_ = facei_;

        // See note in correctAfterParallelTransfer for tetPti_ addressing.
        tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

        // Now the particle is on the receiving side

        // Transform the properties
        if (!receiveCpp.parallel())
        {
            const tensor& T =
            (
                receiveCpp.forwardT().size() == 1
              ? receiveCpp.forwardT()[0]
              : receiveCpp.forwardT()[patchFacei]
            );
            transformProperties(T);
        }
        else if (receiveCpp.separated())
        {
            const vector& s =
            (
                (receiveCpp.separation().size() == 1)
              ? receiveCpp.separation()[0]
              : receiveCpp.separation()[patchFacei]
            );
            transformProperties(-s);
        }
    }
}
//XXXX


template<class TrackData>
void Foam::particle::hitWallPatch
(
    const wallPolyPatch& pp,
    TrackData& td,
    const tetIndices& ti
)
{
    DebugVar(pp.name());
}


template<class TrackData>
void Foam::particle::hitPatch(const polyPatch&, TrackData&)
{}


// ************************************************************************* //
