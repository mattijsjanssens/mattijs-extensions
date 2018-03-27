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

Class
    Foam::fv::skewCorrectedGrad

Description
    Iterative face-based gradient; face values using skew correction.

SourceFiles
    skewCorrectedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef skewCorrectedGrad_H
#define skewCorrectedGrad_H

#include "gradScheme.H"
#include "surfaceInterpolationScheme.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class skewCorrected Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewCorrectedGrad
:
    public fv::gradScheme<Type>
{
    // Private data

        label nIter_;

        //tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        skewCorrectedGrad(const skewCorrectedGrad&);

        //- Disallow default bitwise assignment
        void operator=(const skewCorrectedGrad&);


public:

    //- Runtime type information
    TypeName("skewCorrectedGrad");


    // Constructors

        //- Construct from mesh
        skewCorrectedGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            nIter_(1)
        {
            DebugVar("single iteration");
        }

        //- Construct from mesh and Istream
        skewCorrectedGrad(const fvMesh& mesh, Istream& is)
        :
            gradScheme<Type>(mesh),
            nIter_(is.eof() ? 1 : readLabel(is))
        {
            DebugVar(nIter_);
        }


    // Member Functions

        //- Return the gradient of the given field
        //  calculated using Gauss' theorem on the given surface field
        static
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > gradf
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            const word& name
        );

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;

        //- Correct the boundary values of the gradient using the patchField
        // snGrad functions
        static void correctBoundaryConditions
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "skewCorrectedGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
