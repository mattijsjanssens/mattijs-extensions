/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SubField.H"
#include "transformField.H"

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::patchNeighbourField
(
    const UList<Type>& iF
) const
{
    tmp<Field<Type>> tresult(new Field<Type>(neighbSize()));
    Field<Type>& result = tresult.ref();

    label n = 0;

    for (const label index : AMIIndices_)
    {
        const auto& nbr = neighbPatch(index);
        const labelUList& nbrCells = nbr.faceCells();
        for (const auto celli : nbrCells)
        {
            result[n++] = iF[celli];
        }
    }
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolateUntransformed
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    //if (owner())
    //{
    //    return AMI().interpolateToSource(fld, defaultValues);
    //}
    //else
    //{
    //    return neighbPatch().AMI().interpolateToTarget(fld, defaultValues);
    //}
    tmp<Field<Type>> tresult(new Field<Type>(this->size(), Zero));
    Field<Type>& result = tresult.ref();

    label n = 0;

    for (const label index : AMIIndices_)
    {
        const auto& nbr = neighbPatch(index);
        const auto myAMI(AMI(index));

        if (myAMI.valid())
        {
            // Owner. Interpolate nbr data to here and accumulate
            myAMI().interpolateToSource
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
            const auto nbrAMI(nbr.AMI(neighbIndex(index)));
            if (nbrAMI.valid())
            {
                // Interpolate nbr data to here and accumulate
                nbrAMI().interpolateToTarget
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
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
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

    if (pTraits<Type>::rank == 0)
    {
        return interpolateUntransformed(fld, defaultValues);
    }
    else
    {
        autoPtr<coordSystem::cylindrical> cs(cylindricalCS());
        if (!cs.valid())
        {
            return interpolateUntransformed(fld, defaultValues);
        }
        else
        {
            // Collect all individually transformed nbr fields

            auto tlocalFld(tmp<Field<Type>>::New(fld.size()));
            Field<Type>& localFld = tlocalFld.ref();

            label n = 0;

            for (const label index : AMIIndices_)
            {
                const auto& nbr = neighbPatch(index);

                if (debug)
                {
                    Pout<< "cyclicAMIPolyPatch::interpolate :"
                        << " patch:" << this->name()
                        << " size:" << this->size()
                        << " nbrPatch:" << nbr.name()
                        << " size:" << nbr.size()
                        << endl;
                }

                const SubField<Type> nbrSubFld(fld, nbr.size(), n);
                const Field<Type>& nbrFld(nbrSubFld);

                // Transform to cylindrical coords
                tmp<tensorField> nbrT(cs().R(nbr.faceCentres()));
                const auto tnbrCylFld(Foam::invTransform(nbrT(), nbrFld));
                SubField<Type>(localFld,  nbr.size(), n) = tnbrCylFld;
                n += nbr.size();
            }

            tmp<tensorField> T(cs().R(this->faceCentres()));

            List<Type> localDeflt(defaultValues.size());
            if (defaultValues.size() == size())
            {
                // We get in UList (why? Copied from cyclicAMI). Convert to
                // Field so we can use transformField routines.
                const SubField<Type> defaultSubFld(defaultValues);
                const Field<Type>& defaultFld(defaultSubFld);
                localDeflt = Foam::invTransform(T, defaultFld);
            }

            // Do the actual interpolation and interpolate back to cartesian
            // coords
            return Foam::transform
            (
                T,
                interpolateUntransformed(localFld, localDeflt)
            );
        }
    }
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


template<class Type, class CombineOp>
void Foam::cyclicAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    //- Commented out for now since called with non-primitives (e.g. wallPoint
    //  from FaceCellWave) - these are missing the pTraits<Type>::rank and
    //  Foam::transform
    /*
    autoPtr<coordSystem::cylindrical> cs(cylindricalCS());

    if (cs.valid() && pTraits<Type>::rank > 0)
    {
        const cyclicAMIPolyPatch& nbrPp = this->neighbPatch();

        tmp<tensorField> nbrT(cs().R(nbrPp.faceCentres()));

        result = Foam::invTransform(nbrT, result);
        List<Type> localDeflt(defaultValues.size());
        if (defaultValues.size() == nbrT().size())
        {
            // We get in UList (why? Copied from cyclicAMI). Convert to
            // Field so we can use transformField routines.
            const SubField<Type> defaultSubFld(defaultValues);
            const Field<Type>& defaultFld(defaultSubFld);
            localDeflt = Foam::invTransform(nbrT, defaultFld);
        }

        // Do actual AMI interpolation
        if (owner())
        {
            AMI().interpolateToSource
            (
                fld,
                cop,
                result,
                localDeflt
            );
        }
        else
        {
            neighbPatch().AMI().interpolateToTarget
            (
                fld,
                cop,
                result,
                localDeflt
            );
        }

        // Transform back. Result is now at *this
        const vectorField::subField fc(this->faceCentres());
        result = Foam::transform(cs().R(fc), result);
    }
    else
    */
    {
        if (owner(0))
        {
            AMI(0)().interpolateToSource
            (
                fld,
                cop,
                result,
                defaultValues
            );
        }
        else
        {
            neighbPatch(0).AMI(0)().interpolateToTarget
            (
                fld,
                cop,
                result,
                defaultValues
            );
        }
    }
}


// ************************************************************************* //
