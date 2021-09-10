/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "transform.H"

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

bool Foam::cyclicAMIFvPatch::coupled() const
{
    //return Pstream::parRun() || (this->size() && neighbFvPatch().size());

    if (Pstream::parRun())
    {
        return true;
    }
    else if (!this->size())
    {
        return false;
    }
    else
    {
        return neighbSize() > 0;
    }
}


void Foam::cyclicAMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const scalarField deltas(nf() & coupledFvPatch::delta());

        tmp<scalarField> tnbrDeltas;
        {
            // Collect neighbour deltas

            const auto& pp = cyclicAMIPatch();

            tmp<scalarField> tnd(new scalarField(pp.neighbSize()));
            {
                const labelList& nbrIds = neighbPatchIDs();

                label n = 0;

                scalarField& nd = tnd.ref();
                forAll(nbrIds, index)
                {
                    if (pp.validAMI(index))
                    {
                        const auto& nbrPatch = neighbFvPatch(index);
                        SubField<scalar>(nd, nbrPatch.size(), n) =
                            (nbrPatch.nf() & nbrPatch.coupledFvPatch::delta());
                        n += nbrPatch.size();
                    }
                }
            }

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

        forAll(deltas, facei)
        {
            scalar di = deltas[facei];
            scalar dni = nbrDeltas[facei];

            w[facei] = dni/(di + dni);
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


Foam::tmp<Foam::vectorField> Foam::cyclicAMIFvPatch::delta() const
{
    if (coupled())
    {
        const vectorField patchD(coupledFvPatch::delta());

        tmp<vectorField> tnbrPatchD;
        {
            // Collect neighbour deltas

            const auto& pp = cyclicAMIPatch();

            tmp<vectorField> tpd(new vectorField(pp.neighbSize()));
            {
                const labelList& nbrIds = neighbPatchIDs();

                label n = 0;

                vectorField& pd = tpd.ref();
                forAll(nbrIds, index)
                {
                    if (pp.validAMI(index))
                    {
                        const auto& nbrPatch = neighbFvPatch(index);
                        SubField<vector>(pd, nbrPatch.size(), n) =
                            nbrPatch.coupledFvPatch::delta();
                        n += nbrPatch.size();
                    }
                }
            }

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

        tmp<vectorField> tpdv(new vectorField(patchD.size()));
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


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    tmp<labelField> tfld(new labelField(cyclicAMIPatch().neighbSize()));
    labelField& fld = tfld.ref();

    // Return internal field (e.g. cell agglomeration) in nbr patch index
    const labelList& nbrIds = cyclicAMIPatch().neighbPatchIDs();

    label n = 0;
    forAll(nbrIds, nbri)
    {
        if (cyclicAMIPatch().validAMI(nbri))
        {
            const cyclicAMIFvPatch& nbr = neighbFvPatch(nbri);
            SubField<label>(fld, nbr.size(), n) = nbr.patchInternalField(iF);
            n += nbr.size();
        }
    }

    return tfld;
}


// ************************************************************************* //
