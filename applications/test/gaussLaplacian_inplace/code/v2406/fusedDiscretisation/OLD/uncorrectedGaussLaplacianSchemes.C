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

#include "uncorrectedGaussLaplacianScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvLaplacianScheme(uncorrectedGaussLaplacianScheme)

#define declareFvmLaplacianScalarGamma(Type)                                   \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::fvMatrix<Foam::Type>>                                          \
Foam::fv::uncorrectedGaussLaplacianScheme<Foam::Type, Foam::scalar>::          \
fvmLaplacian                                                                   \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,           \
    const GeometricField<Type, fvPatchField, volMesh>& vf                      \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
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
        if (mesh.fluxRequired(vf.name()))                                      \
        {                                                                      \
            fvm.faceFluxCorrectionPtr() = std::make_unique                     \
            <                                                                  \
                GeometricField<Type, fvsPatchField, surfaceMesh>               \
            >                                                                  \
            (                                                                  \
                gammaMagSf*this->tsnGradScheme_().correction(vf)               \
            );                                                                 \
                                                                               \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    *fvm.faceFluxCorrectionPtr()                               \
                )().primitiveField();                                          \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    gammaMagSf*this->tsnGradScheme_().correction(vf)           \
                )().primitiveField();                                          \
        }                                                                      \
    }                                                                          \
                                                                               \
    return tfvm;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField, Foam::volMesh>> \
Foam::fv::uncorrectedGaussLaplacianScheme<Foam::Type, Foam::scalar>::          \
fvcLaplacian                                                                   \
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
    return tLaplacian;                                                         
}

/*                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField, Foam::volMesh>> \
Foam::fv::uncorrectedGaussLaplacianScheme<Foam::Type, Foam::scalar>::          \
fvcLaplacian                                                                   \
(                                                                              \
    const GeometricField<scalar, fvPatchField, volMesh>& vgamma,               \
    const GeometricField<Type, fvPatchField, volMesh>& vf                      \
)                                                                              \
{                                                                              \
    /* const fvMesh& mesh = this->mesh(); */                                   \
                                                                               \
    auto tgamma(this->tinterpGammaScheme_().interpolate(vgamma));         \
    const auto& gamma = tgamma();                                         \
    return fvcLaplacian(gamma, vf);                                       \
    /* return fvcLaplacian2(vgamma, vf); */ \
}

*/

declareFvmLaplacianScalarGamma(scalar);
//declareFvmLaplacianScalarGamma(vector);
//declareFvmLaplacianScalarGamma(sphericalTensor);
//declareFvmLaplacianScalarGamma(symmTensor);
//declareFvmLaplacianScalarGamma(tensor);


// ************************************************************************* //
