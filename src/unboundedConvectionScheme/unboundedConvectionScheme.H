/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

Class
    Foam::fv::unboundedConvectionScheme

Description
    Wrapper around convectionScheme to bypass the 'bounded' warning.

\*---------------------------------------------------------------------------*/

#ifndef unboundedConvectionScheme_H
#define unboundedConvectionScheme_H

#include "convectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class unboundedConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unboundedConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

        tmp<fv::convectionScheme<Type>> scheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        unboundedConvectionScheme(const unboundedConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const unboundedConvectionScheme&);


public:

    //- Runtime type information
    TypeName("unbounded");


    // Constructors

        //- Construct from flux and Istream
        unboundedConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            scheme_
            (
                fv::convectionScheme<Type>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const surfaceScalarField& phi,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return scheme_().interpolate(phi, vf);
        }

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField& faceFlux,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return scheme_().flux(faceFlux, vf);
        }

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField& faceFlux,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return scheme_().fvmDiv(faceFlux, vf);
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField& faceFlux,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return scheme_().fvcDiv(faceFlux, vf);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
