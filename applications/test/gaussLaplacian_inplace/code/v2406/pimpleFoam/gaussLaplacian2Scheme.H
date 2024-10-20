/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::fv::GaussLaplacianScheme

Group
    grpFvLaplacianSchemes

Description
    Variant of GaussLaplacian that avoids intermediate fields

SourceFiles
    gaussLaplacian2Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef gaussLaplacian2Scheme_H
#define gaussLaplacian2Scheme_H

#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class gaussLaplacian2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class gaussLaplacian2Scheme
:
    public fv::laplacianScheme<Type, GType>
{
    // Private Member Functions

        //- See gaussLaplacianScheme
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- No copy construct
        gaussLaplacian2Scheme(const gaussLaplacian2Scheme&)
            = delete;

        //- No copy assignment
        void operator=(const gaussLaplacian2Scheme&) = delete;


public:

    //- Runtime type information
    TypeName("Gauss2");


    // Constructors

        //- Construct null
        gaussLaplacian2Scheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        gaussLaplacian2Scheme(const fvMesh& mesh, Istream& is)
        :
            laplacianScheme<Type, GType>(mesh, is)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        gaussLaplacian2Scheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}


    //- Destructor
    virtual ~gaussLaplacian2Scheme() = default;


    // Member Functions

        static tmp<fvMatrix<Type>> fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        ////- Helper: grad on component. Move to GaussGrad? Template on
        ////  - cell-access operator?
        ////  - cell-to-face operator?
        ////  - face-to-cell operator?
        //static void gradComponent
        //(
        //    const surfaceScalarField& weights,
        //    const GeometricField<Type, fvPatchField, volMesh>& vf,
        //    const direction cmpt,
        //    GeometricField<Type, fvPatchField, volMesh>& gGrad
        //);

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


    // New: override laplacianScheme::fvcLaplacian for scalar gamma and
    //      scalar,vector volField

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// Use macros to emulate partial-specialisation of the Laplacian functions
// for scalar diffusivity gamma

#define defineFvmLaplacianScalarGamma(Type)                                    \
                                                                               \
template<>                                                                     \
tmp<fvMatrix<Type>> gaussLaplacian2Scheme<Type, scalar>::fvmLaplacian          \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);                                                                             \
                                                                               \
template<>                                                                     \
tmp<GeometricField<Type, fvPatchField, volMesh>>                               \
gaussLaplacian2Scheme<Type, scalar>::fvcLaplacian                              \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);

defineFvmLaplacianScalarGamma(scalar);
defineFvmLaplacianScalarGamma(vector);
defineFvmLaplacianScalarGamma(sphericalTensor);
defineFvmLaplacianScalarGamma(symmTensor);
defineFvmLaplacianScalarGamma(tensor);


template<>
tmp<GeometricField<scalar, fvPatchField, volMesh>>
gaussLaplacian2Scheme<scalar, scalar>::fvcLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>&,
    const GeometricField<scalar, fvPatchField, volMesh>&
);
template<>
tmp<GeometricField<vector, fvPatchField, volMesh>>
gaussLaplacian2Scheme<vector, scalar>::fvcLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>&,
    const GeometricField<vector, fvPatchField, volMesh>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gaussLaplacian2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
