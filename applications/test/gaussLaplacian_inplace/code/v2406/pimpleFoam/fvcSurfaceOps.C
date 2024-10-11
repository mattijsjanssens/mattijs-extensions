/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 M.Janssens
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
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class ResultType, class CellToFaceOp>
void surfaceSum
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& lambdas,
    const CellToFaceOp& cop,
    GeometricField<ResultType, fvPatchField, volMesh>& result,
    const bool doCorrectBoundaryConditions
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // See e.g. surfaceInterpolationScheme<Type>::dotInterpolate

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& lambda = lambdas.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],
                    lambda[facei],
                    vfi[ownCelli],
                    vfi[neiCelli]
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
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
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            pLambda[facei],
                            vfi[pFaceCells[facei]],
                            pnf[facei]
                        )
                    );

                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            scalar(1.0),
                            pvf[facei],
                            pTraits<Type>::one  // not used
                        )
                    );
                    sfi[pFaceCells[facei]] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


template<class Type, class ResultType, class CombineOp>
void GaussOp
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& lambdas,
    const CombineOp& cop,
    GeometricField<ResultType, fvPatchField, volMesh>& result
)
{
    const fvMesh& mesh = vf.mesh();

    // Sum contributions
    surfaceSum(vf, lambdas, cop, result, false);

    auto& sfi = result.primitiveFieldRef();
    sfi /= mesh.V();

    result.correctBoundaryConditions();
}


template<class Type, class ResultType, class CombineOp>
void surfaceOp
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceVectorField& ownLs,
    const surfaceVectorField& neiLs,
    const CombineOp& cop,
    GeometricField<ResultType, fvPatchField, volMesh>& result
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const Type faceVal
            (
                cop
                (
                    Sfi[facei],     // needed?
                    vfi[ownCelli],
                    vfi[neiCelli]
                )
            );
            sfi[ownCelli] += ownLs[facei]*faceVal;
            sfi[neiCelli] -= neiLs[facei]*faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pOwnLs = ownLs.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const Type faceVal
                    (
                        cop
                        (
                            pSf[facei], // needed?
                            vfi[pFaceCells[facei]],
                            pnf[facei]
                        )
                    );

                    sfi[pFaceCells[facei]] += pOwnLs[facei]*faceVal;
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const Type faceVal
                    (
                        cop
                        (
                            pSf[facei],
                            vfi[pFaceCells[facei]],
                            pvf[facei]
                        )
                    );
                    sfi[pFaceCells[facei]] += pOwnLs[facei]*faceVal;
                }
            }
        }
    }

    result.correctBoundaryConditions();
}


template<class Type, class GType, class ResultType, class CellToFaceOp>
void surfaceSnSum
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const surfaceScalarField& gammaWeights,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& deltaCoeffs,
    const CellToFaceOp& cop,
    GeometricField<ResultType, fvPatchField, volMesh>& result,
    const bool doCorrectBoundaryConditions
)
{
    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    const auto& gammai = gamma.primitiveField();
    const auto& vfi = vf.primitiveField();
    auto& sfi = result.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& weights = gammaWeights.primitiveField();
        const auto& dc = deltaCoeffs.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const label ownCelli = P[facei];
            const label neiCelli = N[facei];

            const ResultType faceVal
            (
                cop
                (
                    Sfi[facei],         // area vector

                    weights[facei],     // interpolation weights
                    gammai[ownCelli],
                    gammai[neiCelli],

                    dc[facei],          // delta coefficients
                    vfi[ownCelli],
                    vfi[neiCelli]
                )
            );
            sfi[ownCelli] += faceVal;
            sfi[neiCelli] -= faceVal;
        }
    }


    // Boundary field
    {
        forAll(mesh.boundary(), patchi)
        {
            const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
            const auto& pSf = Sf.boundaryField()[patchi];
            const auto& pvf = vf.boundaryField()[patchi];
            const auto& pdc = deltaCoeffs.boundaryField()[patchi];
            const auto& pweights = gammaWeights.boundaryField()[patchi];
            const auto& pgamma = gamma.boundaryField()[patchi];

            if (pvf.coupled())
            {
                auto tpnf(pvf.patchNeighbourField());
                auto& pnf = tpnf();
                auto tgammanf(pgamma.patchNeighbourField());
                auto& gammanf = tgammanf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Interpolate between owner-side and neighbour-side values
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],         // area vector

                            pweights[facei],
                            gammai[ownCelli],
                            gammanf[facei],

                            pdc[facei],
                            vfi[ownCelli],
                            pnf[facei]
                        )
                    );

                    sfi[ownCelli] += faceVal;
                }
            }
            else
            {
                auto tpnf(pvf.snGrad());
                auto& pnf = tpnf();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    const label ownCelli = pFaceCells[facei];

                    // Use patch value only
                    const ResultType faceVal
                    (
                        cop
                        (
                            pSf[facei],             // area vector

                            scalar(1.0),
                            pgamma[facei],
                            pTraits<GType>::zero,   // not used


                            scalar(1.0),            // use 100% of pnf
                            pTraits<Type>::zero,
                            pnf[facei]
                        )
                    );
                    sfi[ownCelli] += faceVal;
                }
            }
        }
    }

    if (doCorrectBoundaryConditions)
    {
        result.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
