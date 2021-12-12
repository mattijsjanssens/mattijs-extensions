/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "turbulentTemperatureCyclicAMIFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::
turbulentTemperatureCyclicAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    cyclicAMIFvPatchField<scalar>(p, iF),
    temperatureCoupledBase
    (
        patch(),
        "undefined",
        "undefined",
        "undefined-K",
        "undefined-alpha"
    ),
    mappedPatchBase(patch().patch()),
    TnbrName_("undefined-Tnbr")
{}


Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::
turbulentTemperatureCyclicAMIFvPatchScalarField
(
    const turbulentTemperatureCyclicAMIFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMIFvPatchField<scalar>(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mappedPatchBase(patch().patch()),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    thicknessLayer_(ptf.thicknessLayer_.clone(p.patch())),
    kappaLayers_(ptf.kappaLayers_),
    kappaLayer_(ptf.kappaLayer_.clone(p.patch()))
{}


Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::
turbulentTemperatureCyclicAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMIFvPatchField<scalar>(p, iF, dict),
    temperatureCoupledBase(patch(), dict),
    mappedPatchBase(patch().patch(), dict),
    TnbrName_(dict.get<word>("Tnbr"))
{
    // Read list of layers
    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);
    }
    // Read single additional PatchFunction1
    thicknessLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "thicknessLayer",
        dict
    );
    kappaLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "kappaLayer",
        dict
    );
}


Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::
turbulentTemperatureCyclicAMIFvPatchScalarField
(
    const turbulentTemperatureCyclicAMIFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    cyclicAMIFvPatchField<scalar>(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchBase(patch().patch(), wtcsf),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    thicknessLayer_(wtcsf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(wtcsf.kappaLayers_),
    kappaLayer_(wtcsf.kappaLayer_.clone(patch().patch()))
{}


Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::
turbulentTemperatureCyclicAMIFvPatchScalarField
(
    const turbulentTemperatureCyclicAMIFvPatchScalarField& wtcsf
)
:
    cyclicAMIFvPatchField<scalar>(wtcsf),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchBase(patch().patch(), wtcsf),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    thicknessLayer_(wtcsf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(wtcsf.kappaLayers_),
    kappaLayer_(wtcsf.kappaLayer_.clone(patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    cyclicAMIFvPatchField<scalar>::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
    if (thicknessLayer_)
    {
        thicknessLayer_().autoMap(mapper);
        kappaLayer_().autoMap(mapper);
    }
}


void Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    cyclicAMIFvPatchField<scalar>::rmap(ptf, addr);

    const turbulentTemperatureCyclicAMIFvPatchScalarField& tiptf =
        refCast
        <
            const turbulentTemperatureCyclicAMIFvPatchScalarField
        >(ptf);

    temperatureCoupledBase::rmap(tiptf, addr);
    if (thicknessLayer_)
    {
        thicknessLayer_().rmap(tiptf.thicknessLayer_(), addr);
        kappaLayer_().rmap(tiptf.kappaLayer_(), addr);
    }
}


Foam::tmp<Foam::scalarField>
Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::patchDeltaCoeffs() const
{
    // Use patch-normal delta for all non-coupled BCs. TBD: cache
    const vectorField nHat(patch().nf());
    const vectorField delta(nHat*(nHat & (patch().Cf() - patch().Cn())));
    return 1.0/mag(delta);
}


Foam::tmp<Foam::scalarField>
Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::patchKappa
(
    const scalarField& deltaCoeffs,
    const scalarField& Tp
) const
{
    // Get kappa from relevant thermo
    tmp<scalarField> tk(temperatureCoupledBase::kappa(Tp));

    // Optionally modify with explicit resistance
    if (thicknessLayer_ || thicknessLayers_.size())
    {
        scalarField KDelta(tk*deltaCoeffs);

        // Harmonic averaging of kappa*deltaCoeffs
        {
            KDelta = 1.0/KDelta;
            if (thicknessLayer_)
            {
                const scalar t = db().time().timeOutputValue();
                KDelta +=
                    thicknessLayer_().value(t)
                   /kappaLayer_().value(t);
            }
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    KDelta += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
                }
            }
            KDelta = 1.0/KDelta;
        }

        // Update kappa from KDelta
        tk = KDelta/deltaCoeffs;
    }

    return tk;
}


void Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//    // Since we're inside initEvaluate/evaluate there might be processor
//    // comms underway. Change the tag we use.
//    const int oldTag = UPstream::msgType();
//    UPstream::msgType() = oldTag+1;
//
//    // Get the coupling information from the mappedPatchBase
//    const mappedPatchBase& mpp =
//        mappedPatchFieldBase<scalar>::mapper
//        (
//            patch(),
//            this->internalField()
//        );
//
//
//    const scalarField myDeltaCoeffs(patchDeltaCoeffs());
//
//    const scalarField& Tp = *this;
//    const scalarField kappaTp(patchKappa(myDeltaCoeffs, Tp));
//    const tmp<scalarField> myKDelta = kappaTp*myDeltaCoeffs;
//
//
//    scalarField nbrIntFld;
//    scalarField nbrKDelta;
//    if (mpp.sameWorld())
//    {
//        // Same world so lookup
//        const auto& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
//        const label nbrPatchID = mpp.samplePolyPatch().index();
//        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];
//
//        const turbulentTemperatureCyclicAMIFvPatchScalarField&
//        nbrField =
//        refCast
//        <
//        const turbulentTemperatureCyclicAMIFvPatchScalarField
//        >
//        (
//            nbrPatch.lookupPatchField<volScalarField, scalar>
//            (
//                TnbrName_
//            )
//        );
//
//        // Swap to obtain full local values of neighbour K*delta
//        nbrIntFld = nbrField.patchInternalField();
//        nbrKDelta = nbrField.patchKappa(nbrField)*nbrPatch.patchDeltaCoeffs();
//    }
//    else
//    {
//        // Different world so use my region,patch. Distribution below will
//        // do the reordering.
//        nbrIntFld = patchInternalField();
//        nbrKDelta = myKDelta.ref();
//    }
//    distribute(this->internalField().name() + "_value", nbrIntFld);
//    distribute(this->internalField().name() + "_weights", nbrKDelta);
//
//
//    // Both sides agree on
//    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
//    // - gradient    : (temperature-fld)*delta
//    // We've got a degree of freedom in how to implement this in a mixed bc.
//    // (what gradient, what fixedValue and mixing coefficient)
//    // Two reasonable choices:
//    // 1. specify above temperature on one side (preferentially the high side)
//    //    and above gradient on the other. So this will switch between pure
//    //    fixedvalue and pure fixedgradient
//    // 2. specify gradient and temperature such that the equations are the
//    //    same on both sides. This leads to the choice of
//    //    - refGradient = zero gradient
//    //    - refValue = neighbour value
//    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())
//
//    this->refValue() = nbrIntFld;
//    this->refGrad() = 0.0;
//    this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta());
//
//    cyclicAMIFvPatchField<scalar>::updateCoeffs();
//
//    if (debug)
//    {
//        scalar Q = gSum(kappaTp*patch().magSf()*snGrad());
//
//        Info<< patch().boundaryMesh().mesh().name() << ':'
//            << patch().name() << ':'
//            << this->internalField().name() << " <- "
//            << mpp.sampleRegion() << ':'
//            << mpp.samplePatch() << ':'
//            << this->internalField().name() << " :"
//            << " heat transfer rate:" << Q
//            << " walltemperature "
//            << " min:" << gMin(*this)
//            << " max:" << gMax(*this)
//            << " avg:" << gAverage(*this)
//            << endl;
//    }
//
//    // Restore tag
//    UPstream::msgType() = oldTag;
}


//void Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::manipulateMatrix
//(
//    fvMatrix<scalar>& m,
//    const label iMatrix,
//    const direction cmpt
//)
//{
//    FatalErrorInFunction
//        << "This BC does not support energy coupling "
//        << "Use compressible::turbulentTemperatureRadCoupledMixed "
//        << "which has more functionalities and it can handle "
//        << "the assemble coupled option for energy. "
//        << abort(FatalError);
//}
//
//
//Foam::tmp<Foam::Field<Foam::scalar>>
//turbulentTemperatureCyclicAMIFvPatchScalarField::coeffs
//(
//    fvMatrix<scalar>& matrix,
//    const Field<scalar>& coeffs,
//    const label mat
//) const
//{
//    FatalErrorInFunction
//        << "This BC does not support energy coupling "
//        << "Use compressible::turbulentTemperatureRadCoupledMixed "
//        << "which has more functionalities and it can handle "
//        << "the assemble coupled option for energy. "
//        << abort(FatalError);
//    /*
//    const label index(this->patch().index());
//
//    const label nSubFaces(matrix.lduMesh().cellBoundMap()[mat][index].size());
//
//    Field<scalar> mapCoeffs(nSubFaces, Zero);
//
//    label subFaceI = 0;
//    forAll(*this, faceI)
//    {
//
//    }
//    */
//    return tmp<Field<scalar>>(new Field<scalar>());
//}


void Foam::turbulentTemperatureCyclicAMIFvPatchScalarField::write
(
    Ostream& os
) const
{
    cyclicAMIFvPatchField<scalar>::write(os);
    temperatureCoupledBase::write(os);
    mappedPatchBase::write(os);
    os.writeEntry("Tnbr", TnbrName_);
    if (thicknessLayer_)
    {
        thicknessLayer_().writeData(os);
        kappaLayer_().writeData(os);
    }
    if (thicknessLayers_.size())
    {
        thicknessLayers_.writeEntry("thicknessLayers", os);
        kappaLayers_.writeEntry("kappaLayers", os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        turbulentTemperatureCyclicAMIFvPatchScalarField
    );
}

// ************************************************************************* //
