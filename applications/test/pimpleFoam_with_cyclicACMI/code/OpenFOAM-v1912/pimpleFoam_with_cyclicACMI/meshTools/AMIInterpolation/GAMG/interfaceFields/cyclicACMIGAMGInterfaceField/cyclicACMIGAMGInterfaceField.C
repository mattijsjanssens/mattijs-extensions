/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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

#include "cyclicACMIGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"
#include "SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicACMIGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicACMIGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    const cyclicAMILduInterfaceField& p =
        refCast<const cyclicAMILduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(doTransform),
    rank_(rank)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicACMIGAMGInterfaceField::~cyclicACMIGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicACMIGAMGInterfaceField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // Collect all neighbour data (for pairs with valid AMI). See
    // also cyclicAMIPolyPatch::patchNeighbourField)

    solveScalarField pnf(cyclicACMIInterface_.neighbSize());

    label n = 0;
    for (const label index : cyclicACMIInterface_.AMIIndices())
    {
        const auto& nbr = cyclicACMIInterface_.neighbPatch(index);
        const labelUList& nbrCells = nbr.faceCells();
        for (const auto celli : nbrCells)
        {
            pnf[n++] = psiInternal[celli];
        }
    }

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    tmp<solveScalarField> tfld(new solveScalarField(this->size(), Zero));
    solveScalarField& fld = tfld.ref();

    n = 0;
    for (const label index : cyclicACMIInterface_.AMIIndices())
    {
        const auto& nbr = cyclicACMIInterface_.neighbPatch(index);
        const auto myAMI(cyclicACMIInterface_.AMI(index));
        const auto nbrAMI(nbr.AMI(cyclicACMIInterface_.neighbIndex(index)));

        if (myAMI.valid())
        {
            // Owner. Interpolate nbr data to here and accumulate
            myAMI().interpolateToSource
            (
                SubField<solveScalar>(pnf, nbr.size(), n),
                multiplyWeightedOp<solveScalar, plusEqOp<solveScalar>>
                (
                    plusEqOp<solveScalar>()
                ),
                fld
            );
            n += nbr.size();
        }
        else if (nbrAMI.valid())
        {
            // Interpolate nbr data to here and accumulate
            nbrAMI().interpolateToTarget
            (
                SubField<solveScalar>(pnf, nbr.size(), n),
                multiplyWeightedOp<solveScalar, plusEqOp<solveScalar>>
                (
                    plusEqOp<solveScalar>()
                ),
                fld
            );
            n += nbr.size();
        }
    }

    this->addToInternalField(result, !add, coeffs, fld);
}


// ************************************************************************* //
