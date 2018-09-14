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

#include "cyclicAMIFvPatchField.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::cyclicAMIFvPatchField<Foam::scalar>::patchNeighbourField() const
{
    const Field<scalar>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

    Field<scalar> pnf(iField, nbrFaceCells);

    tmp<Field<scalar>> tpnf;
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
        typedef GeometricField<scalar, fvPatchField, volMesh> FieldType;

        const FieldType& fld = dynamic_cast<const FieldType&>
        (
            this->internalField()
        );

const Switch oldCorrected(corrected_);
corrected_ = false;

        typedef typename outerProduct<vector, scalar>::type GradType;
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
            tpnf = cyclicAMIPatch_.interpolate(pnf, gradPnf, UList<scalar>());
        }

corrected_ = oldCorrected;

    }

    if (doTransform())
    {
        tpnf.ref() = transform(forwardT(), tpnf());
    }

    return tpnf;
}


template<>
void Foam::cyclicAMIFvPatchField<Foam::scalar>::updateInterfaceMatrix
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
        typedef GeometricField<scalar, fvPatchField, volMesh> FieldType;

        const FieldType& fld = dynamic_cast<const FieldType&>
        (
            this->internalField()
        );

        const Switch oldCorrected(corrected_);
        corrected_ = false;

        typedef typename outerProduct<vector, scalar>::type GradType;
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
            pnf = cyclicAMIPatch_.interpolate(pnf, gradPnf, UList<scalar>());
        }
    }

    // Multiply the field by coefficients and add into the result
    const labelUList& faceCells = cyclicAMIPatch_.faceCells();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// template<>
// void Foam::cyclicAMIFvPatchField<Foam::scalar>::updateInterfaceMatrix
// (
//     Field<scalar>& result,
//     const Field<scalar>& psiInternal,
//     const scalarField& coeffs,
//     const Pstream::commsTypes
// ) const
// {
//     const labelUList& nbrFaceCells =
//         cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();
// 
//     Field<scalar> pnf(psiInternal, nbrFaceCells);
// 
//     // Transform according to the transformation tensors
//     transformCoupleField(pnf);
// 
//     if (!corrected_)
//     {
//         if (cyclicAMIPatch_.applyLowWeightCorrection())
//         {
//             Field<scalar> pif(psiInternal, cyclicAMIPatch_.faceCells());
//             pnf = cyclicAMIPatch_.interpolate(pnf, pif);
//         }
//         else
//         {
//             pnf = cyclicAMIPatch_.interpolate(pnf);
//         }
//     }
//     else
//     {
//         typedef GeometricField<scalar, fvPatchField, volMesh> FieldType;
// 
//         const FieldType& fld = dynamic_cast<const FieldType&>
//         (
//             this->internalField()
//         );
// 
//         const Switch oldCorrected(corrected_);
//         corrected_ = false;
// 
//         typedef typename outerProduct<vector, scalar>::type GradType;
//         tmp<GeometricField<GradType, fvPatchField, volMesh>> tgradFld
//         (
//             fvc::grad(fld)
//         );
//         Field<GradType> gradPnf(tgradFld(), nbrFaceCells);
// 
//         corrected_ = oldCorrected;
// 
//         if (cyclicAMIPatch_.applyLowWeightCorrection())
//         {
//             Field<scalar> pif(psiInternal, cyclicAMIPatch_.faceCells());
//             pnf = cyclicAMIPatch_.interpolate(pnf, gradPnf, pif);
//         }
//         else
//         {
//             pnf = cyclicAMIPatch_.interpolate(pnf, gradPnf);
//         }
//     }
// 
//     // Multiply the field by coefficients and add into the result
//     const labelUList& faceCells = cyclicAMIPatch_.faceCells();
// 
//     forAll(faceCells, elemI)
//     {
//         result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
//     }
// }


template<>
void Foam::cyclicAMIFvPatchField<Foam::scalar>::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
    os.writeKeyword("corrected") << corrected_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
