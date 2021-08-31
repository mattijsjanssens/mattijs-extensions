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
    const labelList& nbrPatchIds = neighbPatchIDs();

    tmp<Field<Type>> tresult(new Field<Type>(neighbSize()));
    Field<Type>& result = tresult.ref();

    label n = 0;

    forAll(nbrPatchIds, index)
    {
        if (validAMI(index))
        {
            const auto& nbrPp = neighbPatch(index);
            const labelUList& nbrCells = nbrPp.faceCells();
            for (const auto celli : nbrCells)
            {
                result[n++] = iF[celli];
            }
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
    const labelList& nbrPatchIds = neighbPatchIDs();

    // fld = nbr values
    // defaultValues = local fall-back values

    if
    (
        debug
     && (
            fld.size() != neighbSize()
         || (defaultValues.size() && defaultValues.size() != size())
        )
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

    label n = 0;

    forAll(nbrPatchIds, index)
    {
        const auto& nbr = neighbPatch(index);

        if (AMI(index).valid())
        {
            //Pout<< "** owner:" << this->name() << " size:" << this->size()
            //    << " interpolateToSource from " << nbr.name()
            //    << " size:" << nbr.size()
            //    << " collated field:" << fld.size() << " starting at:" << n
            //    << " into " << result.size() << endl;

            // Owner. Interpolate nbr data to here and accumulate
            AMI(index)().interpolateToSource
            (
                SubField<Type>(fld, nbr.size(), n),
                multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
                result,
                defaultValues
            );
            n += nbr.size();
        }
        else
        {
            // Check if slave has valid AMI
            const label myIndex = nbr.neighbPatchIDs().find(this->index());
            if (nbr.AMI(myIndex).valid())
            {
                //Pout<< "** neighour:" << nbr.name()
                //    << " size:" << nbr.size()
                //    << " interpolateToSource from " << nbr.name()
                //    << " size:" << nbr.size()
                //    << " collated field:" << fld.size()
                //    << " starting at:" << n
                //    << " into " << result.size() << endl;

                // Interpolate nbr data to here and accumulate
                nbr.AMI(myIndex)().interpolateToTarget
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


// ************************************************************************* //
