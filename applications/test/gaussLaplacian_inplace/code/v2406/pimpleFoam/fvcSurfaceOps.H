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

InNamespace
    Foam::fvc

Description
    Surface integrate surfaceField creating a volField.
    Surface sum a surfaceField creating a volField.

SourceFiles
    fvcSurfaceIntegrate.C

\*---------------------------------------------------------------------------*/


#ifndef fvcSurfaceIntegrate_H
#define fvcSurfaceIntegrate_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    //- Interpolate to face (using cop) and apply Gauss. Note: uses V(),
    //  not Vsc()
    template<class Type, class ResultType, class CombineOp>
    void GaussOp
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const surfaceScalarField& lambdas,
        const CombineOp& cop,
        GeometricField<ResultType, fvPatchField, volMesh>& result
    );

    //- Interpolate to face (using cop) and apply vectors
    template<class Type, class ResultType, class CombineOp>
    void surfaceOp
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const surfaceVectorField& ownLs,
        const surfaceVectorField& neiLs,
        const CombineOp& cop,
        GeometricField<ResultType, fvPatchField, volMesh>& result
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcSurfaceOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
