/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "skewCorrectedGrad.H"
#include "skewCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::skewCorrectedGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        fv::gaussGrad<Type>::gradf(ssf, name)
    );

    return tgGrad;
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::skewCorrectedGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

DebugVar(name);
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf
    (
        //tinterpScheme_().interpolate(vsf)
        linearInterpolate(vsf)
    );
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf();


    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        gradf(ssf, name)
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    const skewCorrectionVectors& skv = skewCorrectionVectors::New(vsf.mesh());
    for (label i = 0; i < nIter_; i++)
    {
DebugVar(i);
        tmp<GeometricField<GradType, fvsPatchField, surfaceMesh>> tsgGrad
        (
            linearInterpolate(gGrad)
        );

DebugVar(i);
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tcorr
        (
            skv()
          & tsgGrad
        );
        tcorr.ref().dimensions().reset(vsf.dimensions());
DebugVar(i);

        const GeometricField<GradType, fvPatchField, volMesh> oldGrad(gGrad);

DebugVar(i);

        gGrad =
            0.9*gGrad
          + 0.1*fv::gaussGrad<Type>::gradf(tcorr+ssf, name);

DebugVar(i);
        const GeometricField<GradType, fvPatchField, volMesh> resid
        (
            gGrad-oldGrad
        );
        Pout<< "ITer:" << i << " residual:" << gAverage(resid)
            << " max:" << gMax(resid) << endl;
    }

    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
void Foam::fv::skewCorrectedGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
DebugVar(vsf.name());


    typename GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >::Boundary& gGradbf = gGrad.boundaryFieldRef();

    forAll(vsf.boundaryField(), patchi)
    {
        if (!vsf.boundaryField()[patchi].coupled())
        {
            const vectorField n
            (
                vsf.mesh().Sf().boundaryField()[patchi]
              / vsf.mesh().magSf().boundaryField()[patchi]
            );

            gGradbf[patchi] += n *
            (
                vsf.boundaryField()[patchi].snGrad()
              - (n & gGradbf[patchi])
            );
        }
     }
}


// ************************************************************************* //
