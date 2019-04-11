/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "relaxedNonOrthoGaussLaplacianScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvLaplacianScheme(relaxedNonOrthoGaussLaplacianScheme)

#define declareFvmLaplacianScalarGamma(Type)                                   \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::fvMatrix<Foam::Type>>                                          \
Foam::fv::relaxedNonOrthoGaussLaplacianScheme<Foam::Type, Foam::scalar>::      \
fvmLaplacian                                                                   \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,           \
    const GeometricField<Type, fvPatchField, volMesh>& vf                      \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SType;            \
                                                                               \
    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf              \
    (                                                                          \
        gamma*mesh.magSf()                                                     \
    );                                                                         \
                                                                               \
    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected                         \
    (                                                                          \
        gammaMagSf,                                                            \
        this->tsnGradScheme_().deltaCoeffs(vf),                                \
        vf                                                                     \
    );                                                                         \
    fvMatrix<Type>& fvm = tfvm.ref();                                          \
                                                                               \
    if (this->tsnGradScheme_().corrected())                                    \
    {                                                                          \
        tmp<SType> tfaceFluxCorrection                                         \
        (                                                                      \
            gammaMagSf                                                         \
           *this->tsnGradScheme_().correction(vf)                              \
        );                                                                     \
                                                                               \
        tmp<SType> trelaxedCorrection(new SType(tfaceFluxCorrection()));       \
                                                                               \
        if (oldCorrection_.valid())                                            \
        {                                                                      \
/* Pout<< "oldCorr:" << gAverage(oldCorrection_()) << " newCorr:" << gAverage(tfaceFluxCorrection()) << endl; */ \
            Pout<< "Underrelaxing " << vf.name()                               \
                << " non-ortho by 0.5" << endl;                                \
            trelaxedCorrection.ref() *= 0.5;                                   \
            trelaxedCorrection.ref() += 0.5*oldCorrection_();                  \
        }                                                                      \
                                                                               \
        tmp<Field<Type>> tcorr                                                 \
        (                                                                      \
            mesh.V()                                                           \
           *fvc::div                                                           \
            (                                                                  \
                trelaxedCorrection()                                           \
            )().primitiveField()                                               \
        );                                                                     \
                                                                               \
        fvm.source() -= tcorr();                                               \
                                                                               \
Pout<< "Transferring non-ortho " << vf.name() << endl;                         \
        oldCorrection_ = tfaceFluxCorrection;                                  \
                                                                               \
        if (mesh.fluxRequired(vf.name()))                                      \
        {                                                                      \
            fvm.faceFluxCorrectionPtr() = trelaxedCorrection.ptr();            \
        }                                                                      \
    }                                                                          \
                                                                               \
    return tfvm;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField, Foam::volMesh>> \
Foam::fv::relaxedNonOrthoGaussLaplacianScheme<Foam::Type, Foam::scalar>::fvcLaplacian         \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,           \
    const GeometricField<Type, fvPatchField, volMesh>& vf                      \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian                \
    (                                                                          \
        fvc::div(gamma*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())         \
    );                                                                         \
                                                                               \
    tLaplacian.ref().rename                                                    \
    (                                                                          \
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'                    \
    );                                                                         \
                                                                               \
    return tLaplacian;                                                         \
}


declareFvmLaplacianScalarGamma(scalar);
declareFvmLaplacianScalarGamma(vector);
declareFvmLaplacianScalarGamma(sphericalTensor);
declareFvmLaplacianScalarGamma(symmTensor);
declareFvmLaplacianScalarGamma(tensor);


// ************************************************************************* //
