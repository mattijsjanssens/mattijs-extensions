/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "cyclicAMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "transform.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicAMIFvPatch, polyPatch);
    addNamedToRunTimeSelectionTable
    (
        fvPatch,
        cyclicAMIFvPatch,
        polyPatch,
        cyclicPeriodicAMI
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::cyclicAMIFvPatch::newInternalProcFaces
// (
//     label& newFaces,
//     label& newProcFaces
// ) const
// {
//     const labelListList& addSourceFaces = AMI().srcAddress();
//
//     // Add new faces as many weights for AMI
//     forAll (addSourceFaces, faceI)
//     {
//         const labelList& nbrFaceIs = addSourceFaces[faceI];
//
//         forAll (nbrFaceIs, j)
//         {
//             label nbrFaceI = nbrFaceIs[j];
//
//             if (nbrFaceI < neighbPatch().size())
//             {
//                 // local faces
//                 newFaces++;
//             }
//             else
//             {
//                 // Proc faces
//                 newProcFaces++;
//             }
//         }
//     }
// }


bool Foam::cyclicAMIFvPatch::coupled() const
{
    return
        Pstream::parRun()
     || !this->boundaryMesh().mesh().time().processorCase();
}


void Foam::cyclicAMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const scalarField deltas(nf() & coupledFvPatch::delta());

Pout<< "** cyclicAMIFvPatch::makeWeights : patch:" << this->name()
    << " deltas:" << deltas << endl;
Pout<< "    neighbSize:" << neighbSize() << " AMIIndices:" << AMIIndices()
    << endl;

        tmp<scalarField> tnbrDeltas;
        {
            // Collect neighbour deltas

            tmp<scalarField> tnd(new scalarField(neighbSize()));
            {
                //const labelList& nbrIds = neighbPatchIDs();

                label n = 0;

                scalarField& nd = tnd.ref();
                for (const label index : AMIIndices())
                {
                    const auto& nbr = neighbFvPatch(index);

Pout<< "nbr:" << nbr.name() << " nf:" << nbr.nf()
    << " delta:" << nbr.coupledFvPatch::delta() << endl;

                    SubField<scalar>(nd, nbr.size(), n) =
                        (nbr.nf() & nbr.coupledFvPatch::delta());
                    n += nbr.size();
                }
            }

        Pout<< "** cyclicAMIFvPatch::makeWeights : patch:" << this->name()
            << " NON-INTERPLATED nbrDeltas:"
            << tnd() << endl;

            // Interpolate to *this
            if (applyLowWeightCorrection())
            {
                tnbrDeltas = interpolate(tnd, scalarField(this->size(), 1.0));
            }
            else
            {
                tnbrDeltas = interpolate(tnd);
            }
        }

        const scalarField& nbrDeltas = tnbrDeltas();

Pout<< "** cyclicAMIFvPatch::makeWeights : patch:" << this->name()
    << " nbrDeltas:" << nbrDeltas << endl;

        forAll(deltas, facei)
        {
            // Note use of mag
            scalar di = mag(deltas[facei]);
            scalar dni = mag(nbrDeltas[facei]);

            w[facei] = dni/(di + dni);
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


void Foam::cyclicAMIFvPatch::makeDeltaCoeffs(scalarField& coeffs) const
{
    // Apply correction to default coeffs
}


void Foam::cyclicAMIFvPatch::makeNonOrthoDeltaCoeffs(scalarField& coeffs) const
{
    // Apply correction to default coeffs
    //coeffs = Zero;
}


void Foam::cyclicAMIFvPatch::makeNonOrthoCorrVectors(vectorField& vecs) const
{
    // Apply correction to default vectors
    //vecs = Zero;
}


Foam::tmp<Foam::vectorField> Foam::cyclicAMIFvPatch::delta() const
{
    if (coupled())
    {
        const vectorField patchD(coupledFvPatch::delta());

        tmp<vectorField> tnbrPatchD;
        //if (applyLowWeightCorrection())
        //{
        //    tnbrPatchD =
        //        interpolate
        //        (
        //            nbrPatch.coupledFvPatch::delta(),
        //            vectorField(this->size(), Zero)
        //        );
        //}
        //else
        //{
        //    tnbrPatchD = interpolate(nbrPatch.coupledFvPatch::delta());
        //}
        {
            // Collect neighbour deltas

            tmp<vectorField> tpd(new vectorField(neighbSize()));
            {
                //const labelList& nbrIds = neighbPatchIDs();

                label n = 0;

                vectorField& pd = tpd.ref();
                for (const label index : AMIIndices())
                {
                    const auto& nbr = neighbFvPatch(index);
                    SubField<vector>(pd, nbr.size(), n) =
                        nbr.coupledFvPatch::delta();
                    n += nbr.size();
                }
            }

Pout<< "** cyclicAMIFvPatch::delta : NON-INTEROLATED patch:" << this->name()
    << " tpd:" << tpd() << endl;

            // Interpolate to *this
            if (applyLowWeightCorrection())
            {
                tnbrPatchD = interpolate(tpd, vectorField(this->size(), Zero));
            }
            else
            {
                tnbrPatchD = interpolate(tpd);
            }
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


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    return patchInternalField(internalData, faceCells);
}


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::internalFieldTransfer
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
        const auto& nbr = neighbFvPatch(index);
        SubField<label>(fld, nbr.size(), n) = nbr.patchInternalField(iF);
        n += nbr.size();
    }

    return tfld;
}


void Foam::cyclicAMIFvPatch::movePoints()
{
    if (!owner(0) || !cyclicAMIPolyPatch_.createAMIFaces())
    {
        // Only manipulating patch face areas and mesh motion flux if the AMI
        // creates additional faces
        return;
    }

    // Update face data based on values set by the AMI manipulations
    const_cast<vectorField&>(Sf()) = cyclicAMIPolyPatch_.faceAreas();
    const_cast<vectorField&>(Cf()) = cyclicAMIPolyPatch_.faceCentres();
    const_cast<scalarField&>(magSf()) = mag(Sf());

    const cyclicAMIFvPatch& nbr = neighbPatch(0);
    const_cast<vectorField&>(nbr.Sf()) = nbr.cyclicAMIPatch().faceAreas();
    const_cast<vectorField&>(nbr.Cf()) = nbr.cyclicAMIPatch().faceCentres();
    const_cast<scalarField&>(nbr.magSf()) = mag(nbr.Sf());


    // Set consitent mesh motion flux
    // TODO: currently maps src mesh flux to tgt - update to
    // src = src + mapped(tgt) and tgt = tgt + mapped(src)?

    const fvMesh& mesh = boundaryMesh().mesh();
    surfaceScalarField& meshPhi = const_cast<fvMesh&>(mesh).setPhi();
    surfaceScalarField::Boundary& meshPhiBf = meshPhi.boundaryFieldRef();

    if (cyclicAMIPolyPatch_.owner(0))
    {
        scalarField& phip = meshPhiBf[patch().index()];
        forAll(phip, facei)
        {
            const face& f = cyclicAMIPolyPatch_.localFaces()[facei];

            // Note: using raw point locations to calculate the geometric
            // area - faces areas are currently scaled by the AMI weights
            // (decoupled from mesh points)
            const scalar geomArea = f.mag(cyclicAMIPolyPatch_.localPoints());

            const scalar scaledArea = magSf()[facei];
            phip[facei] *= scaledArea/geomArea;
        }

        scalarField srcMeshPhi(phip);
        if (AMI(0)().distributed())
        {
            AMI(0)().srcMap().distribute(srcMeshPhi);
        }

        const labelListList& tgtToSrcAddr = AMI(0)().tgtAddress();
        scalarField& nbrPhip = meshPhiBf[nbr.index()];

        forAll(tgtToSrcAddr, tgtFacei)
        {
            // Note: now have 1-to-1 mapping so tgtToSrcAddr[tgtFacei] is size 1
            const label srcFacei = tgtToSrcAddr[tgtFacei][0];
            nbrPhip[tgtFacei] = -srcMeshPhi[srcFacei];
        }

        DebugInfo
            << "patch:" << patch().name()
            << " sum(area):" << gSum(magSf())
            << " min(mag(faceAreas):" << gMin(magSf())
            << " sum(meshPhi):" << gSum(phip) << nl
            << " sum(nbrMeshPhi):" << gSum(nbrPhip) << nl
            << endl;
    }
}

// ************************************************************************* //
