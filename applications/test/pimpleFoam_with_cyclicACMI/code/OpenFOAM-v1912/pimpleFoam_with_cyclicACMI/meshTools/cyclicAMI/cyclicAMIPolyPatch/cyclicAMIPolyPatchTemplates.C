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

#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::patchNeighbourField
(
    const UList<Type>& iF
) const
{
    const labelList& nbrIds = neighbPatchIDs();

    label n = 0;
    for (const label nbrId : nbrIds)
    {
        n += this->boundaryMesh()[nbrId].size();
    }

    tmp<Field<Type>> tresult(new Field<Type>(n));
    Field<Type>& result = tresult.ref();

    n = 0;
    for (const label nbrId : nbrIds)
    {
        const labelUList& nbrCells = this->boundaryMesh()[nbrId].faceCells();
        for (const auto celli : nbrCells)
        {
            result[n++] = iF[celli];
        }
    }
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
//    const cyclicAMIPolyPatch& ownPatch
//    (
//        owner()
//      ? *this
//      : neighbPatch()
//    );
//
//    const labelList& nbrIds = ownPatch.neighbPatchIDs();
    const labelList& nbrIds = neighbPatchIDs();

    // fld = nbr values
    // defaultValues = local fall-back values

    if
    (
        fld.size() != neighbSize()
     || (defaultValues.size() && defaultValues.size() != size())
    )
    {
        FatalErrorInFunction << "Field size:" << fld.size()
            << " neighbour patches:" << neighbPatchIDs()
            << " total size:" << neighbSize()
            << " defaultValues:" << defaultValues.size()
            << " local size:" << size()
            << exit(FatalError);
    }

    tmp<Field<Type>> tresult(new Field<Type>(this->size(), Zero));
    Field<Type>& result = tresult.ref();

    if (owner())
    {
        //return AMI().interpolateToSource(fld, defaultValues);

        label n = 0;
        forAll(nbrIds, nbri)
        {
            const cyclicAMIPolyPatch& nbr = neighbPatch(nbri);

            // Interpolate nbr data to here and accumulate
            AMI(nbri).interpolateToSource
            (
                SubField<Type>(fld, nbr.size(), n),
                multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
                result,
                defaultValues
            );
            n += nbr.size();
        }
    }
    else
    {
        //return neighbPatch().AMI().interpolateToTarget(fld, defaultValues);

        label n = 0;
        forAll(nbrIds, nbri)
        {
            const cyclicAMIPolyPatch& nbr = neighbPatch(nbri);

            if (nbr.owner())
            {
                // Find the AMI that points to me
                label myId = nbr.neighbPatchIDs().find(this->index());

                // Interpolate nbr data to here and accumulate
                nbr.AMI(myId).interpolateToTarget
                (
                    SubField<Type>(fld, nbr.size(), n),
                    multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
                    result,
                    defaultValues
                );
                n += nbr.size();
            }
        }
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), defaultValues);
}


//template<class Type, class CombineOp>
//void Foam::cyclicAMIPolyPatch::interpolate
//(
//    const UList<Type>& fld,
//    const CombineOp& cop,
//    List<Type>& result,
//    const UList<Type>& defaultValues
//) const
//{
//    if (owner())
//    {
//        AMI().interpolateToSource
//        (
//            fld,
//            cop,
//            result,
//            defaultValues
//        );
//    }
//    else
//    {
//        neighbPatch().AMI().interpolateToTarget
//        (
//            fld,
//            cop,
//            result,
//            defaultValues
//        );
//    }
//}


// ************************************************************************* //
