/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "mappedPatchFieldBase.H"
#include "mappedPatchBase.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Type Foam::mappedPatchFieldBase<Type>::getAverage
(
    const dictionary& dict,
    const bool mandatory
)
{
    if (mandatory)
    {
        return dict.get<Type>("average");
    }

    return Zero;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const word& fieldName,
    const bool setAverage,
    const Type average,
    const word& interpolationScheme
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(fieldName),
    setAverage_(setAverage),
    average_(average),
    interpolationScheme_(interpolationScheme)
{}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const dictionary& dict
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_
    (
        dict.template lookupOrDefault<word>
        (
            "field",
            patchField_.internalField().name()
        )
    ),
    setAverage_(dict.lookupOrDefault("setAverage", false)),
    average_(getAverage(dict, setAverage_)),
    interpolationScheme_(interpolationCell<Type>::typeName)
{
    if (mapper_.mode() == mappedPatchBase::NEARESTCELL)
    {
        dict.readEntry("interpolationScheme", interpolationScheme_);
    }
}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(patchField_.internalField().name()),
    setAverage_(false),
    average_(Zero),
    interpolationScheme_(interpolationCell<Type>::typeName)
{}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchFieldBase<Type>& mapper
)
:
    mapper_(mapper.mapper_),
    patchField_(mapper.patchField_),
    fieldName_(mapper.fieldName_),
    setAverage_(mapper.setAverage_),
    average_(mapper.average_),
    interpolationScheme_(mapper.interpolationScheme_)
{}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const mappedPatchFieldBase<Type>& base
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(base.fieldName_),
    setAverage_(base.setAverage_),
    average_(base.average_),
    interpolationScheme_(base.interpolationScheme_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&
Foam::mappedPatchFieldBase<Type>::sampleField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (mapper_.sameRegion())
    {
        if (fieldName_ == patchField_.internalField().name())
        {
            // Optimisation: bypass field lookup
            return
                dynamic_cast<const fieldType&>
                (
                    patchField_.internalField()
                );
        }
        else
        {
            const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
            return thisMesh.template lookupObject<fieldType>(fieldName_);
        }
    }

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());
    return nbrMesh.template lookupObject<fieldType>(fieldName_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchFieldBase<Type>::mappedField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    // Result of obtaining remote values
    auto tnewValues = tmp<Field<Type>>::New();
    auto& newValues = tnewValues.ref();

    switch (mapper_.mode())
    {
        case mappedPatchBase::NEARESTCELL:
        {
            const mapDistribute& distMap = mapper_.map();

            if (interpolationScheme_ != interpolationCell<Type>::typeName)
            {
                // Send back sample points to the processor that holds the cell
                vectorField samples(mapper_.samplePoints());
                distMap.reverseDistribute
                (
                    (
                        mapper_.sameRegion()
                      ? thisMesh.nCells()
                      : nbrMesh.nCells()
                    ),
                    point::max,
                    samples
                );

                auto interpolator =
                    interpolation<Type>::New
                    (
                        interpolationScheme_,
                        sampleField()
                    );

                const auto& interp = *interpolator;

                newValues.setSize(samples.size(), pTraits<Type>::max);
                forAll(samples, celli)
                {
                    if (samples[celli] != point::max)
                    {
                        newValues[celli] = interp.interpolate
                        (
                            samples[celli],
                            celli
                        );
                    }
                }
            }
            else
            {
                newValues = sampleField();
            }

            distMap.distribute(newValues);

            break;
        }
        case mappedPatchBase::NEARESTPATCHFACE:
        case mappedPatchBase::NEARESTPATCHFACEAMI:
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

            if (nbrPatchID < 0)
            {
                FatalErrorInFunction
                 << "Unable to find sample patch " << mapper_.samplePatch()
                 << " in region " << mapper_.sampleRegion()
                 << " for patch " << patchField_.patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues = nbrField.boundaryField()[nbrPatchID];
            mapper_.distribute(newValues);

            break;
        }
        case mappedPatchBase::NEARESTFACE:
        {
            Field<Type> allValues(nbrMesh.nFaces(), Zero);

            const fieldType& nbrField = sampleField();

            for (const fvPatchField<Type>& pf : nbrField.boundaryField())
            {
                label faceStart = pf.patch().start();

                forAll(pf, facei)
                {
                    allValues[faceStart++] = pf[facei];
                }
            }

            mapper_.distribute(allValues);
            newValues.transfer(allValues);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown sampling mode: " << mapper_.mode() << nl
                << abort(FatalError);
        }
    }

    if (setAverage_)
    {
        Type averagePsi =
            gSum(patchField_.patch().magSf()*newValues)
           /gSum(patchField_.patch().magSf());

        if (mag(averagePsi)/mag(average_) > 0.5)
        {
            newValues *= mag(average_)/mag(averagePsi);
        }
        else
        {
            newValues += (average_ - averagePsi);
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    return tnewValues;
}


template<class Type>
void Foam::mappedPatchFieldBase<Type>::write(Ostream& os) const
{
    os.writeEntry("field", fieldName_);

    if (setAverage_)
    {
        os.writeEntry("setAverage", "true");
        os.writeEntry("average", average_);
    }

    os.writeEntry("interpolationScheme", interpolationScheme_);
}


// ************************************************************************* //
