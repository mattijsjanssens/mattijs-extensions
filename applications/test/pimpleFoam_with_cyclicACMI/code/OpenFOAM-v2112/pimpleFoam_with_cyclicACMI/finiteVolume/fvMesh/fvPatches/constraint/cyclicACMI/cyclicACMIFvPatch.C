/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "cyclicACMIFvPatch.H"
#include "fvMesh.H"
#include "transform.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicACMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::cyclicACMIFvPatch::updateAreas() const
{
    // Give AMI chance to update itself
    bool updated = cyclicACMIPolyPatch_.updateAreas();

    //if (!cyclicACMIPolyPatch_.owner())
    //{
    //    return updated;
    //}

    if (updated || !cyclicACMIPolyPatch_.upToDate(areaTime_))
    {
        if (debug)
        {
            Pout<< "cyclicACMIFvPatch::updateAreas() : updating fv areas for "
                << name() << " and " << this->nonOverlapPatch().name()
                << endl;
        }

        const fvPatch& nonOverlapPatch = this->nonOverlapPatch();
        //const cyclicACMIFvPatch& nbrACMI = neighbPatch();
        //const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapPatch();

        resetPatchAreas(*this);
        resetPatchAreas(nonOverlapPatch);
        //resetPatchAreas(nbrACMI);
        //resetPatchAreas(nbrNonOverlapPatch);

        updated = true;

        // Mark my data to be up to date with ACMI polyPatch level
        cyclicACMIPolyPatch_.setUpToDate(areaTime_);
    }
    return updated;
}


void Foam::cyclicACMIFvPatch::resetPatchAreas(const fvPatch& fvp) const
{
    const_cast<vectorField&>(fvp.Sf()) = fvp.patch().faceAreas();
    const_cast<vectorField&>(fvp.Cf()) = fvp.patch().faceCentres();
    const_cast<scalarField&>(fvp.magSf()) = mag(fvp.patch().faceAreas());

    DebugPout
        << fvp.patch().name() << " area:" << sum(fvp.magSf()) << endl;
}


void Foam::cyclicACMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const scalarField deltas(nf() & coupledFvPatch::delta());

        //const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();
        //
        //// These deltas are of the cyclic part alone - they are
        //// not affected by the amount of overlap with the nonOverlapPatch
        //scalarField nbrDeltas
        //(
        //    interpolate
        //    (
        //        nbrPatch.nf() & nbrPatch.coupledFvPatch::delta()
        //    )
        //);

        // Collect nbr delta
        scalarField nbd(neighbSize());
        {
            label n = 0;
            for (const label index : AMIIndices())
            {
                //const auto& nbr = neighbFvPatch(index);
                const auto& nbr = neighbPatch(index);
                SubField<scalar>(nbd, nbr.size(), n) =
                    (nbr.nf() & nbr.coupledFvPatch::delta());
                n += nbr.size();
            }
        }

        // These deltas are of the cyclic part alone - they are
        // not affected by the amount of overlap with the nonOverlapPatch
        scalarField nbrDeltas(interpolate(nbd));


        const scalar tol = cyclicACMIPolyPatch::tolerance();


        forAll(deltas, facei)
        {
            scalar di = mag(deltas[facei]);
            scalar dni = mag(nbrDeltas[facei]);

            if (dni < tol)
            {
                // Avoid zero weights on disconnected faces. This value
                // will be weighted with the (zero) face area so will not
                // influence calculations.
                w[facei] = 1.0;
            }
            else
            {
                w[facei] = dni/(di + dni);
            }
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicACMIFvPatch::cyclicACMIFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    coupledFvPatch(patch, bm),
    cyclicACMILduInterface(),
    cyclicACMIPolyPatch_(refCast<const cyclicACMIPolyPatch>(patch)),
    areaTime_
    (
        IOobject
        (
            "areaTime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    )
{
    areaTime_.eventNo() = -1;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::cyclicACMIFvPatch::newInternalProcFaces
// (
//     label& newFaces,
//     label& newProcFaces
// ) const
// {
//     const List<labelList>& addSourceFaces =
//         cyclicACMIPolyPatch_.AMI().srcAddress();
//
//     const scalarField& fMask = cyclicACMIPolyPatch_.srcMask();
//
//     // Add new faces as many weights for AMI
//     forAll (addSourceFaces, faceI)
//     {
//         if (fMask[faceI] > cyclicACMIPolyPatch_.tolerance_)
//         {
//             const labelList& nbrFaceIs = addSourceFaces[faceI];
//
//             forAll (nbrFaceIs, j)
//             {
//                 label nbrFaceI = nbrFaceIs[j];
//
//                 if (nbrFaceI < neighbPatch().size())
//                 {
//                     // local faces
//                     newFaces++;
//                 }
//                 else
//                 {
//                     // Proc faces
//                     newProcFaces++;
//                 }
//             }
//         }
//     }
// }


// Foam::refPtr<Foam::labelListList>
// Foam::cyclicACMIFvPatch::mapCollocatedFaces() const
// {
//     const scalarField& fMask = cyclicACMIPolyPatch_.srcMask();
//     const labelListList& srcFaces = cyclicACMIPolyPatch_.AMI().srcAddress();
//     labelListList dOverFaces;
//
//     dOverFaces.setSize(srcFaces.size());
//     forAll (dOverFaces, faceI)
//     {
//         if (fMask[faceI] > cyclicACMIPolyPatch_.tolerance_)
//         {
//             dOverFaces[faceI].setSize(srcFaces[faceI].size());
//
//             forAll (dOverFaces[faceI], subFaceI)
//             {
//                 dOverFaces[faceI][subFaceI] = srcFaces[faceI][subFaceI];
//             }
//         }
//     }
//     return refPtr<labelListList>(new labelListList(dOverFaces));
// }


bool Foam::cyclicACMIFvPatch::coupled() const
{
    return Pstream::parRun() || (this->size() && this->neighbSize());
}


Foam::tmp<Foam::vectorField> Foam::cyclicACMIFvPatch::delta() const
{
    if (coupled())
    {
        const vectorField patchD(coupledFvPatch::delta());

        //const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();
        //vectorField nbrPatchD(interpolate(nbrPatch.coupledFvPatch::delta()));

        tmp<vectorField> tnbrPatchD;
        {
            // Collect neighbour deltas

            vectorField pd(neighbSize());
            {
                label n = 0;
                for (const label index : AMIIndices())
                {
                    //const auto& nbr = neighbFvPatch(index);
                    const auto& nbr = neighbPatch(index);

                    SubField<vector>(pd, nbr.size(), n) =
                        nbr.coupledFvPatch::delta();
                    n += nbr.size();
                }
            }
            tnbrPatchD = interpolate(pd);
        }
        const vectorField& nbrPatchD = tnbrPatchD();

        auto tpdv = tmp<vectorField>::New(patchD.size());
        vectorField& pdv = tpdv.ref();

        // do the transformation if necessary
        if (parallel())
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - dni;
            }
        }
        else
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - transform(forwardT()[0], dni);
            }
        }

        return tpdv;
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    return patchInternalField(internalData, faceCells);
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    //return neighbFvPatch().patchInternalField(iF);

    tmp<labelField> tfld(new labelField(neighbSize()));
    labelField& fld = tfld.ref();

    // Return internal field (e.g. cell agglomeration) in nbr patch index
    label n = 0;
    for (const label index : AMIIndices())
    {
        //const auto& nbr = neighbFvPatch(index);
        const auto& nbr = neighbPatch(index);
        SubField<label>(fld, nbr.size(), n) = nbr.patchInternalField(iF);
        n += nbr.size();
    }
    return tfld;
}


//XXXX
void Foam::cyclicACMIFvPatch::resetPatchPhi
(
    const cyclicACMIFvPatch& fvp,
    const scalarField& mask,
    const labelListList& addr
)
{
    const fvMesh& mesh = boundaryMesh().mesh();
    const auto& points = mesh.points();
    auto& meshPhi = const_cast<fvMesh&>(mesh).setPhi();
    auto& meshPhiBf = meshPhi.boundaryFieldRef();
    const auto& pp = fvp.cyclicACMIPatch();

    scalarField& phip = meshPhiBf[pp.index()];
    forAll(phip, facei)
    {
        if (addr[facei].empty())
        {
            phip[facei] = 0.0;
        }
        else
        {
            const face& fAMI = pp[facei];

            // Note: using raw point locations to calculate the geometric
            // area - faces areas are currently scaled (decoupled from
            // mesh points)
            const scalar geomArea = fAMI.mag(points);
            phip[facei] *= fvp.magSf()[facei]/geomArea;
        }
    }

    const auto& nonOverlapPatch = fvp.nonOverlapPatch();
    auto& phiNonOverlapp = meshPhiBf[nonOverlapPatch.patch().index()];
    forAll(phiNonOverlapp, facei)
    {
        const scalar w = 1.0 - mask[facei];
        phiNonOverlapp[facei] *= w;
    }
}
//XXXX

void Foam::cyclicACMIFvPatch::movePoints()
{
//    if (!cyclicACMIPolyPatch_.owner())
//    {
//        return;
//    }


    if (!cyclicACMIPolyPatch_.upToDate(areaTime_))
    {
        if (debug)
        {
            Pout<< "cyclicACMIFvPatch::movePoints() : updating fv areas for "
                << name() << " and " << this->nonOverlapPatch().name()
                << endl;
        }


        // Set the patch face areas to be consistent with the changes made
        // at the polyPatch level


        //const fvMesh& mesh = boundaryMesh().mesh();
        //surfaceScalarField& meshPhi = const_cast<fvMesh&>(mesh).setPhi();
        //surfaceScalarField::Boundary& meshPhiBf = meshPhi.boundaryFieldRef();


        // Make sure fv areas are consistent with poly level
        resetPatchAreas(*this);
        resetPatchAreas(this->nonOverlapPatch());

        for (const label index : AMIIndices())
        {
            const auto& nbrACMI = neighbPatch(index);
            const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapPatch();

            // Make sure fv areas are consistent with poly level
            resetPatchAreas(nbrACMI);
            resetPatchAreas(nbrNonOverlapPatch);
        }

        // Scale meshphi using adapted fv areas
        for (const label index : AMIIndices())
        {
            const auto& nbrACMI = neighbPatch(index);

            // Scale the mesh flux

            const auto myAMI(AMI(index));

            //const fvPatch& nonOverlapPatch = this->nonOverlapPatch();
            //const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapPatch();
            //const labelListList& newSrcAddr = AMI()().srcAddress();
            //const labelListList& newTgtAddr = AMI()().tgtAddress();
            //
            //// Note: phip and phiNonOverlap will be different sizes if new faces
            //// have been added
            //scalarField& phip = meshPhiBf[cyclicACMIPolyPatch_.index()];
            //scalarField& phiNonOverlapp =
            //    meshPhiBf[nonOverlapPatch.patch().index()];
            //
            //const auto& points = mesh.points();
            //
            //forAll(phip, facei)
            //{
            //    if (newSrcAddr[facei].empty())
            //    {
            //        // AMI patch with no connection to other coupled faces
            //        phip[facei] = 0.0;
            //    }
            //    else
            //    {
            //        // Scale the mesh flux according to the area fraction
            //        const face& fAMI = cyclicACMIPolyPatch_[facei];
            //
            //        // Note: using raw point locations to calculate the geometric
            //        // area - faces areas are currently scaled (decoupled from
            //        // mesh points)
            //        const scalar geomArea = fAMI.mag(points);
            //        phip[facei] *= magSf()[facei]/geomArea;
            //    }
            //}
            //
            //forAll(phiNonOverlapp, facei)
            //{
            //    const scalar w = 1.0 - cyclicACMIPolyPatch_.mask()[facei];
            //    phiNonOverlapp[facei] *= w;
            //}


            if (myAMI.valid())
            {
                resetPatchPhi
                (
                    *this,
                    cyclicACMIPolyPatch_.mask(),
                    myAMI().srcAddress()
                );
                resetPatchPhi
                (
                    nbrACMI,
                    nbrACMI.cyclicACMIPatch().mask(),
                    myAMI().tgtAddress()
                );
            }
            else
            {
                //const cyclicACMIPolyPatch& nbrPatch = nbrACMI.cyclicACMIPatch();
                //scalarField& nbrPhip = meshPhiBf[nbrPatch.index()];
                //scalarField& nbrPhiNonOverlapp =
                //    meshPhiBf[nbrNonOverlapPatch.patch().index()];
                //
                //forAll(nbrPhip, facei)
                //{
                //    if (newTgtAddr[facei].empty())
                //    {
                //        nbrPhip[facei] = 0.0;
                //    }
                //    else
                //    {
                //        const face& fAMI = nbrPatch[facei];
                //
                //        // Note: using raw point locations to calculate the geometric
                //        // area - faces areas are currently scaled (decoupled from
                //        // mesh points)
                //        const scalar geomArea = fAMI.mag(points);
                //        nbrPhip[facei] *= nbrACMI.magSf()[facei]/geomArea;
                //    }
                //}
                //
                //forAll(nbrPhiNonOverlapp, facei)
                //{
                //    const scalar w = 1.0 - cyclicACMIPolyPatch_.tgtMask()[facei];
                //    nbrPhiNonOverlapp[facei] *= w;
                //}

                const auto nbrAMI(nbrACMI.AMI(neighbIndex(index)));
                resetPatchPhi
                (
                    *this,
                    cyclicACMIPolyPatch_.mask(),
                    nbrAMI().tgtAddress()
                );
                resetPatchPhi
                (
                    nbrACMI,
                    nbrACMI.cyclicACMIPatch().mask(),
                    nbrAMI().srcAddress()
                );
            }
        }

        // Mark my data to be up to date with ACMI polyPatch level
        cyclicACMIPolyPatch_.setUpToDate(areaTime_);
    }
}

// ************************************************************************* //
