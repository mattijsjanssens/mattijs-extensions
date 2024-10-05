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
#include "gaussDivScheme2.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
template<class Type>
void gaussDivScheme2<Type>::fvcDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& lambdas,
    GeometricField
    <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>& fld
)
{
    typedef typename innerProduct<vector, Type>::type DivType;

    const fvMesh& mesh = vf.mesh();
    const auto& Sf = mesh.Sf();
    const auto& P = mesh.owner();
    const auto& N = mesh.neighbour();

    auto& fldi = fld.primitiveFieldRef();

    // Internal field
    {
        const auto& Sfi = Sf.primitiveField();
        const auto& vfi = vf.primitiveField();
        const auto& lambda = lambdas.primitiveField();

        for (label facei=0; facei<P.size(); facei++)
        {
            const DivType Sfssf =
                Sfi[facei]
              & (lambda[facei]*(vfi[P[facei]] - vfi[N[facei]]) + vfi[N[facei]]);

            fldi[P[facei]] += Sfssf;
            fldi[N[facei]] -= Sfssf;
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
                const auto& pnf = tpnf();
                auto tpif(pvf.patchInternalField());
                const auto& pif = tpif();

                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    fldi[pFaceCells[facei]] +=
                        pSf[facei]
                      & lerp(pif[facei], pnf[facei], pLambda[facei]);
                }
            }
            else
            {
                for (label facei=0; facei<pFaceCells.size(); facei++)
                {
                    fldi[pFaceCells[facei]] += pSf[facei] & pvf[facei];
                }
            }
        }
    }

    fldi /= mesh.Vsc();
    fld.correctBoundaryConditions();
}
*/


template<class Type>
tmp
<
    GeometricField
    <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
>
gaussDivScheme2<Type>::fvcDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    typedef GeometricField<DivType, fvPatchField, volMesh> DivFieldType;

    tmp<DivFieldType> tDiv
    (
        new DivFieldType
        (
            IOobject
            (
                "div(" + vf.name() + ')',
                vf.instance(),
                vf.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            dimensioned<DivType>(vf.dimensions()/dimLength, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    DivFieldType& div = tDiv.ref();

    //fvcDiv(vf, this->tinterpScheme_().weights(vf), div);

    const auto interpolator = [&]
    (
        const vector& area,
        const scalar lambda,
        const Type& ownVal,
        const Type& neiVal
    ) -> DivType
    {
        return area & (lambda*(ownVal - neiVal) + neiVal);
    };

    fvc::GaussOp
    (
        vf,
        this->tinterpScheme_().weights(vf),
        interpolator,
        div
    );

    return tDiv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
