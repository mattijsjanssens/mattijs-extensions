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

#include "fvcSurfaceOps.H"
#include "gaussGrad2.H"
#include "uncorrectedGaussLaplacianScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
uncorrectedGaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected
(
    const surfaceScalarField& gammaMagSf,
    const surfaceScalarField& deltaCoeffs,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.upper() = deltaCoeffs.primitiveField()*gammaMagSf.primitiveField();
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[patchi];
        const fvsPatchScalarField& pDeltaCoeffs =
            deltaCoeffs.boundaryField()[patchi];

        if (pvf.coupled())
        {
            fvm.internalCoeffs()[patchi] =
                pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            fvm.boundaryCoeffs()[patchi] =
               -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
        }
        else
        {
            fvm.internalCoeffs()[patchi] = pGamma*pvf.gradientInternalCoeffs();
            fvm.boundaryCoeffs()[patchi] = -pGamma*pvf.gradientBoundaryCoeffs();
        }
    }

    return tfvm;
}


//XXXXXX
template<class Type, class GType>
void uncorrectedGaussLaplacianScheme<Type, GType>::gradComponent
(
    const surfaceScalarField& weights,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const direction cmpt,
    GeometricField<Type, fvPatchField, volMesh>& gGrad
)
{
    gGrad = Zero;

    const auto interpolator = [&]
    (
        const vector& area,
        const scalar lambda,
        const Type& ownVal,
        const Type& neiVal
    ) -> Type
    {
        return area*(lambda*(ownVal[cmpt] - neiVal[cmpt]) + neiVal[cmpt]);
    };

    fvc::GaussOp
    (
        vf,
        weights,
        interpolator,
        gGrad
    );

    gaussGrad2<Type>::correctBoundaryConditions(vf, gGrad);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
uncorrectedGaussLaplacianScheme<Type, GType>::gammaSnGradCorr
(
    const surfaceVectorField& SfGammaCorr,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tgammaSnGradCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "gammaSnGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            SfGammaCorr.dimensions()
           *vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    tgammaSnGradCorr.ref().oriented() = SfGammaCorr.oriented();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tgammaSnGradCorr.ref().replace
        (
            cmpt,
            fvc::dotInterpolate(SfGammaCorr, fvc::grad(vf.component(cmpt)))
        );
    }

    return tgammaSnGradCorr;
}
template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
uncorrectedGaussLaplacianScheme<Type, GType>::gammaSnGradCorr
(
    const surfaceScalarField& weights,
    const surfaceVectorField& SfGammaCorr,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tgammaSnGradCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "gammaSnGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            SfGammaCorr.dimensions()
           *vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    tgammaSnGradCorr.ref().oriented() = SfGammaCorr.oriented();


    GeometricField<Type, fvPatchField, volMesh> gradCmptFld
    (
        IOobject
        (
            vf.name() + ".component()",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        vf.dimensions()
    );

    //const auto tweights(this->tinterpGammaScheme_().weights(gamma));
    //const auto& weights = tweights();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        gradComponent
        (
            weights,
            vf,
            cmpt,
            gradCmptFld
        );

        tgammaSnGradCorr.ref().replace
        (
            cmpt,
            //fvc::dotInterpolate(SfGammaCorr, fvc::grad(vf.component(cmpt)))
            fvc::dotInterpolate(SfGammaCorr, gradCmptFld)
        );
    }

    return tgammaSnGradCorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
uncorrectedGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    tmp<FieldType> tresult
    (
        new FieldType
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                vf.instance(),
                vf.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            dimensioned<Type>(vf.dimensions()/dimArea, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    FieldType& result = tresult.ref();


    // Note: cannot use fvc::GaussOp since specialised handling on boundary.
    // Maybe bypass for processor boundaries?

    const auto tdeltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
    const auto& lambdas = tdeltaCoeffs();


    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& lambda = lambdas.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownFacei = P[facei];
            const label neiFacei = N[facei];
            const Type faceVal
            (
                mag(Sfi[facei])
              * lambda[facei]
              * (vfi[neiFacei]-vfi[ownFacei])
            );
            sfi[ownFacei] += faceVal;
            sfi[neiFacei] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pLambda = lambdas.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.snGrad(pLambda));
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const Type faceVal(mag(pSf[facei])*pnf[facei]);
                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
            else
            {
                auto tpnf(pvf.snGrad());
                auto& pnf = tpnf();
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Use patch value only
                    const Type faceVal(mag(pSf[facei])*pnf[facei]);
                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
        }
    }

    sfi /= mesh.V();
    result.correctBoundaryConditions();

    return tresult;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
uncorrectedGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    tmp<FieldType> tresult
    (
        new FieldType
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                vf.instance(),
                vf.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            dimensioned<Type>(gamma.dimensions()*vf.dimensions()/dimArea, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    FieldType& result = tresult.ref();

    const auto tweights(this->tinterpGammaScheme_().weights(gamma));
    const auto& weights = tweights();
    const auto tdeltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
    const auto& dcs = tdeltaCoeffs();

/*
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& dc = dcs.primitiveField();
        const auto& weightsi = weights.primitiveField();
        const auto& vfi = vf.primitiveField();
        const auto& gammai = gamma.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownFacei = P[facei];
            const label neiFacei = N[facei];

            const auto ownVal = vfi[ownFacei];
            const auto neiVal = vfi[neiFacei];

            const vector Sf = Sfi[facei];
            const scalar magSf(mag(Sf));
            const vector Sn(Sf/magSf);

            const auto ownGamma = gammai[ownFacei];
            const auto neiGamma = gammai[neiFacei];

            // Interpolated value
            const GType faceGamma(weightsi[facei]*(ownGamma-neiGamma)+neiGamma);

            // Normal & tangential component of faceGamma
            const vector SfGamma(Sf & faceGamma);
            const scalar SfGammaSn(SfGamma & Sn);
            //const vector SfGammCorr(SfGamma - SfGammaSn*Sn);

            const Type snGrad(dc[facei] * (neiVal-ownVal));

            // Difference across face multiplied by area magnitude
            const Type faceVal(magSf * SfGammaSn * snGrad);

            sfi[ownFacei] += faceVal;
            sfi[neiFacei] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pdc = dcs.boundaryField()[patchi];
            const auto& pgamma = gamma.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.snGrad(pdc));
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const vector Sf = pSf[facei];
                    const scalar magSf(mag(Sf));
                    const vector Sn(Sf/magSf);

                    const auto faceGamma = pgamma[facei];

                    // Normal & tangential component of faceGamma
                    const vector SfGamma(Sf & faceGamma);
                    const scalar SfGammaSn(SfGamma & Sn);
                    //const vector SfGammCorr(SfGamma - SfGammaSn*Sn);

                    const Type faceVal(mag(pSf[facei])*SfGammaSn*pnf[facei]);
                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
            else
            {
                auto tpnf(pvf.snGrad());
                auto& pnf = tpnf();
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Use patch value only
                    const vector Sf = pSf[facei];
                    const scalar magSf(mag(Sf));
                    const vector Sn(Sf/magSf);

                    const auto faceGamma = pgamma[facei];

                    // Normal & tangential component of faceGamma
                    const vector SfGamma(Sf & faceGamma);
                    const scalar SfGammaSn(SfGamma & Sn);

                    const Type faceVal(mag(pSf[facei])*SfGammaSn*pnf[facei]);
                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
        }
    }
*/

    const auto snGrad = [&]
    (
        const vector& Sf,

        const scalar weight,
        const GType& ownGamma,
        const GType& neiGamma,

        const scalar dc,
        const Type& ownVal,
        const Type& neiVal
    ) -> Type
    {
        const scalar magSf(mag(Sf));
        const vector Sn(Sf/magSf);
        const GType faceGamma(weight*(ownGamma-neiGamma)+neiGamma);

        // Normal & tangential component of faceGamma
        const vector SfGamma(Sf & faceGamma);
        const scalar SfGammaSn(SfGamma & Sn);

        const Type snGrad(dc*(neiVal-ownVal));
        return magSf*SfGammaSn*snGrad;
    };

    fvc::surfaceSnSum
    (
        gamma,
        weights,

        vf,
        dcs,

        snGrad,

        result,
        false
    );

    result.primitiveFieldRef() /= vf.mesh().V();
    result.correctBoundaryConditions();

    return tresult;
}
//XXXXXXXX
//template<class Type, class GType>
//tmp<GeometricField<Type, fvPatchField, volMesh>>
//uncorrectedGaussLaplacianScheme<Type, GType>::fvcLaplacian
//(
//    const GeometricField<scalar, fvPatchField, volMesh>& gamma,
//    const GeometricField<Type, fvPatchField, volMesh>& vf
//)
//{
//    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
//
//    tmp<FieldType> tresult
//    (
//        new FieldType
//        (
//            IOobject
//            (
//                "laplacian(" + vf.name() + ')',
//                vf.instance(),
//                vf.mesh(),
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            vf.mesh(),
//            dimensioned<Type>(vf.dimensions()/dimArea, Zero),
//            fvPatchFieldBase::extrapolatedCalculatedType()
//        )
//    );
//    FieldType& result = tresult.ref();
//
//    const auto tweights(this->tinterpGammaScheme_().weights(gamma));
//    const auto& weights = tweights();
//
//
//    // Note: cannot use fvc::GaussOp since specialised handling on boundary.
//    // Maybe bypass for processor boundaries?
//
//    const auto tdeltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
//    const auto& dcs = tdeltaCoeffs();
//
//
//    const fvMesh& mesh = vf.mesh();
//    const auto& Sf = mesh.Sf();
//    const auto& P = mesh.owner();
//    const auto& N = mesh.neighbour();
//
//    const auto& gammai = gamma.primitiveField();
//    const auto& vfi = vf.primitiveField();
//    auto& sfi = result.primitiveFieldRef();
//
//    // Internal field
//    {
//        const auto& Sfi = Sf.primitiveField();
//        const auto& dc = dcs.primitiveField();
//        const auto& weightsi = weights.primitiveField();
//
//        for (label facei=0; facei<P.size(); facei++)
//        {
//            const label ownFacei = P[facei];
//            const label neiFacei = N[facei];
//
//            const auto ownVal = vfi[ownFacei];
//            const auto neiVal = vfi[neiFacei];
//
//            const vector Sf = Sfi[facei];
//            const scalar magSf(mag(Sf));
//
//            const auto ownGamma = gammai[ownFacei];
//            const auto neiGamma = gammai[neiFacei];
//
//            // Interpolated value
//            const GType faceGamma(weightsi[facei]*(ownGamma-neiGamma)+neiGamma);
//
//            const Type snGrad(dc[facei] * (neiVal-ownVal));
//
//            // Difference across face multiplied by area magnitude
//            const Type faceVal(magSf*faceGamma*snGrad);
//
//            sfi[ownFacei] += faceVal;
//            sfi[neiFacei] -= faceVal;
//        }
//    }
//
//
//    // Boundary field
//    {
//        forAll(mesh.boundary(), patchi)
//        {
//            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
//            const auto& pSf = Sf.boundaryField()[patchi];
//            const auto& pvf = vf.boundaryField()[patchi];
//            const auto& pdc = dcs.boundaryField()[patchi];
//            const auto& pgamma = gamma.boundaryField()[patchi];
//
//            if (pvf.coupled())
//            {
//                auto tpnf(pvf.snGrad(pdc));
//                auto& pnf = tpnf();
//
//                for (label facei=0; facei<pFaceCells.size(); facei++)
//                {
//                    const vector Sf = pSf[facei];
//                    const scalar magSf(mag(Sf));
//
//                    const auto faceGamma = pgamma[facei];
//                    const Type faceVal(magSf*faceGamma*pnf[facei]);
//                    sfi[pFaceCells[facei]] += faceVal;
//                }
//            }
//            else
//            {
//                auto tpnf(pvf.snGrad());
//                auto& pnf = tpnf();
//                for (label facei=0; facei<pFaceCells.size(); facei++)
//                {
//                    // Use patch value only
//                    const vector Sf = pSf[facei];
//                    const scalar magSf(mag(Sf));
//
//                    const auto faceGamma = pgamma[facei];
//                    const Type faceVal(magSf*faceGamma*pnf[facei]);
//                    sfi[pFaceCells[facei]] += faceVal;
//                }
//            }
//        }
//    }
//
//////XXXXXXXXX
////    const auto snGrad = [&]
////    (
////        const vector& Sf,
////
////        const scalar weight,
////        const GType ownGamma,
////        const GType neiGamma,
////
////        const scalar dc,
////        const Type& ownVal,
////        const Type& neiVal
////    ) -> Type
////    {
////        const GType faceGamma(weight*(ownGamma-neiGamma)+neiGamma);
////        const Type snGrad(dc*(neiVal-ownVal));
////        return mag(Sf)*faceGamma*snGrad;
////    };
////
////    fvc::surfaceSnSum
////    (
////        gamma,
////        weights,
////        vf,
////        dcs,
////        snGrad,
////        result,
////        false
////    );
//////XXXXXXXXX
//
//
//    sfi /= mesh.V();
//    result.correctBoundaryConditions();
//
//    return tresult;
//}
//XXXXXXXX

template<class Type, class GType>
tmp<fvMatrix<Type>>
uncorrectedGaussLaplacianScheme<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());

    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
    (
        SfGamma & Sn
    );
    const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);

    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected
    (
        SfGammaSn,
        this->tsnGradScheme_().deltaCoeffs(vf),
        vf
    );
    fvMatrix<Type>& fvm = tfvm.ref();

//    const auto tweights(this->tinterpGammaScheme_().weights(gamma));
//    const auto& weights = tweights();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfaceFluxCorrection
        = gammaSnGradCorr(SfGammaCorr, vf);

    if (this->tsnGradScheme_().corrected())
    {
        tfaceFluxCorrection.ref() +=
            SfGammaSn*this->tsnGradScheme_().correction(vf);
    }

    fvm.source() -= mesh.V()*fvc::div(tfaceFluxCorrection())().primitiveField();

    if (mesh.fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr(tfaceFluxCorrection.ptr());
    }

    return tfvm;
}


// TBD. Same as gaussLaplacianScheme::fvcLaplacian
template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
uncorrectedGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());
    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
    (
        SfGamma & Sn
    );
    const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);

    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian
    (
        fvc::div
        (
            SfGammaSn*this->tsnGradScheme_().snGrad(vf)
          + gammaSnGradCorr(SfGammaCorr, vf)
        )
    );

    tLaplacian.ref().rename
    (
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );

    return tLaplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
