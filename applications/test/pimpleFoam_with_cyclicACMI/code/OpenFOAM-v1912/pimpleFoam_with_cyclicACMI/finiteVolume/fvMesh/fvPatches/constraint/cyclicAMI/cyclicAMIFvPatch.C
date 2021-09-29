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

Pout<< "** cyclicAMIFvPatch::delta : patch:" << this->name()
    << " patchD:" << patchD << endl;

        tmp<vectorField> tnbrPatchD;
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

Pout<< "** cyclicAMIFvPatch::delta : patch:" << this->name()
    << " nbrPatchD:" << nbrPatchD << endl;


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


// ************************************************************************* //
