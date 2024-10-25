/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024 M. Janssens
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

#include "fvcSurfaceOps.H"
#include "gaussConvectionScheme2.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const surfaceInterpolationScheme<Type>&
gaussConvectionScheme2<Type>::interpScheme() const
{
    return tinterpScheme_();
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
gaussConvectionScheme2<Type>::interpolate
(
    const surfaceScalarField&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return tinterpScheme_().interpolate(vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
gaussConvectionScheme2<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
gaussConvectionScheme2<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    //fvm.lower() = -weights.primitiveField()*faceFlux.primitiveField();
    multiplySubtract
    (
        fvm.lower(),
        weights.primitiveField(),
        faceFlux.primitiveField()
    );

    //fvm.upper() = fvm.lower() + faceFlux.primitiveField();
    add(fvm.upper(), fvm.lower(), faceFlux.primitiveField());

    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& psf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

        fvm.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }

    if (tinterpScheme_().corrected())
    {
        fvm += fvc::surfaceIntegrate(faceFlux*tinterpScheme_().correction(vf));
    }

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
gaussConvectionScheme2<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const fvMesh& mesh = vf.mesh();

    tmp<FieldType> tConvection
    (
        new FieldType
        (
            IOobject
            (
                "convection(" + faceFlux.name() + ',' + vf.name() + ')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(vf.dimensions()/dimLength, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );

    if (this->tinterpScheme_().corrected())
    {
        const auto tfaceCorr(this->tinterpScheme_().correction(vf));
        auto& faceCorr = tfaceCorr();

        const auto interpolator = [&]
        (
            const vector& area,
            const scalar lambda,

            const Type& ownVal,
            const Type& neiVal,

            const scalar& faceVal,
            const Type& correction

        ) -> Type
        {
            return faceVal*((lambda*(ownVal - neiVal) + neiVal) + correction);
        };

        fvc::surfaceSum
        (
            // interpolation factors for volume field
            this->tinterpScheme_().weights(vf),

            // volume field(s)
            vf,

            // surface field(s)
            faceFlux,
            faceCorr,

            // operation
            interpolator,

            tConvection.ref(),
            false
        );
    }
    else
    {
        const auto interpolator = [&]
        (
            const vector& area,
            const scalar lambda,
            const Type& ownVal,
            const Type& neiVal,

            const scalar& faceVal
        ) -> Type
        {
            return faceVal*(lambda*(ownVal - neiVal) + neiVal);
        };

        fvc::surfaceSum
        (
            tinterpScheme_().weights(vf),

            vf,

            faceFlux,

            interpolator,
            tConvection.ref(),
            false
        );
    }

    tConvection.ref().primitiveFieldRef() /= mesh.Vsc();

    tConvection.ref().correctBoundaryConditions();

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
