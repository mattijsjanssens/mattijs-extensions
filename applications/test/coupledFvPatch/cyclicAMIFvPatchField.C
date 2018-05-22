/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    cyclicAMILduInterfaceField(),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    corrected_(false)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    cyclicAMILduInterfaceField(),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    corrected_(dict.lookup("corrected"))
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value") && this->coupled())
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMILduInterfaceField(),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    corrected_(ptf.corrected_)
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf
)
:
    coupledFvPatchField<Type>(ptf),
    cyclicAMILduInterfaceField(),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    corrected_(ptf.corrected_)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    cyclicAMILduInterfaceField(),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    corrected_(ptf.corrected_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::coupled() const
{
    return cyclicAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(iField, nbrFaceCells);

    tmp<Field<Type>> tpnf;
    if (!corrected_)
    {
        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            tpnf = cyclicAMIPatch_.interpolate
            (
                pnf,
                this->patchInternalField()()
            );
        }
        else
        {
            tpnf = cyclicAMIPatch_.interpolate(pnf);
        }
    }
    else
    {
        typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

        const FieldType& fld = dynamic_cast<const FieldType&>
        (
            this->internalField()
        );

const Switch oldCorrected(corrected_);
corrected_ = false;

        typedef typename outerProduct<vector, Type>::type GradType;
        tmp<GeometricField<GradType, fvPatchField, volMesh>> tgradFld
        (
            fvc::grad(fld)
        );

//DebugVar(fld);
//DebugVar(tgradFld());


        Field<GradType> gradPnf(tgradFld(), nbrFaceCells);

Pout<< "Field:" << fld.name() << " patch:" << this->patch().name()
    << " size:" << this->patch().size()
    << " nbrFld:" << pnf
    << " correcting with:" << gradPnf << endl;

        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            tpnf = cyclicAMIPatch_.interpolate
            (
                pnf,
                gradPnf,
                this->patchInternalField()()
            );
        }
        else
        {
            tpnf = cyclicAMIPatch_.interpolate(pnf, gradPnf);
        }

corrected_ = oldCorrected;

    }

    if (doTransform())
    {
        tpnf.ref() = transform(forwardT(), tpnf());
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicAMIFvPatchField<Type>&
Foam::cyclicAMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const cyclicAMIFvPatchField<Type>>
    (
        fld.boundaryField()[cyclicAMIPatch_.neighbPatchID()]
    );
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    if (!corrected_)
    {
        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            scalarField pif(psiInternal, cyclicAMIPatch_.faceCells());
            pnf = cyclicAMIPatch_.interpolate(pnf, pif);
        }
        else
        {
            pnf = cyclicAMIPatch_.interpolate(pnf);
        }
    }
    else
    {
        typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

        const FieldType& fld = dynamic_cast<const FieldType&>
        (
            this->internalField()
        );

        const Switch oldCorrected(corrected_);
        corrected_ = false;

        typedef typename outerProduct<vector, Type>::type GradType;
        tmp<GeometricField<GradType, fvPatchField, volMesh>> tgradFld
        (
            fvc::grad(fld)
        );
        Field<GradType> gradPnf(tgradFld(), nbrFaceCells);

        corrected_ = oldCorrected;

        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            scalarField pif(psiInternal, cyclicAMIPatch_.faceCells());
            pnf = cyclicAMIPatch_.interpolate(pnf, gradPnf, pif);
        }
        else
        {
            pnf = cyclicAMIPatch_.interpolate(pnf, gradPnf);
        }
    }

    // Multiply the field by coefficients and add into the result
    const labelUList& faceCells = cyclicAMIPatch_.faceCells();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    if (!corrected_)
    {
        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            Field<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());
            pnf = cyclicAMIPatch_.interpolate(pnf, pif);
        }
        else
        {
            pnf = cyclicAMIPatch_.interpolate(pnf);
        }
    }
    else
    {
        typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

        const FieldType& fld = dynamic_cast<const FieldType&>
        (
            this->internalField()
        );

        const Switch oldCorrected(corrected_);
        corrected_ = false;

        typedef typename outerProduct<vector, Type>::type GradType;
        tmp<GeometricField<GradType, fvPatchField, volMesh>> tgradFld
        (
            fvc::grad(fld)
        );
        Field<GradType> gradPnf(tgradFld(), nbrFaceCells);

        corrected_ = oldCorrected;

        if (cyclicAMIPatch_.applyLowWeightCorrection())
        {
            Field<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());
            pnf = cyclicAMIPatch_.interpolate(pnf, gradPnf, pif);
        }
        else
        {
            pnf = cyclicAMIPatch_.interpolate(pnf, gradPnf);
        }
    }

    // Multiply the field by coefficients and add into the result
    const labelUList& faceCells = cyclicAMIPatch_.faceCells();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
    os.writeKeyword("corrected") << corrected_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
