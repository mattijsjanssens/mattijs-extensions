/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::simpleSolver

Description

SourceFiles
    simpleSolver.C

\*---------------------------------------------------------------------------*/

#ifndef simpleSolver_H
#define simpleSolver_H

#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "IOMRFZoneList.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class meshToMesh;

/*---------------------------------------------------------------------------*\
                         Class simpleSolver Declaration
\*---------------------------------------------------------------------------*/

class simpleSolver
{
    // Private data

        //- Runtime
        const Time& runTime_;

        //- The mesh
        const fvMesh& mesh_;

        //- Controls for SIMPLE
        simpleControl simple_;

        //- Pressure field
        volScalarField p_;

        //- Velocity field
        volVectorField U_;

        //- Flux field
        surfaceScalarField phi_;

        //- Pressure reference cell
        label pRefCell_;

        //- Pressure reference value
        scalar pRefValue_;

        //- Laminar transport model
        singlePhaseTransportModel laminarTransport_;

        //- Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence_;

        //- Multiple reference frames
        IOMRFZoneList MRF_;

        //- FV-options
        fv::options fvOptions_;

        //- Previous velocity field
        tmp<volVectorField> U0_;

        //- Previous pressure field
        tmp<volScalarField> p0_;


        tmp<volScalarField> k0_;
        tmp<volScalarField> omega0_;
        tmp<volScalarField> nut0_;

        //- Velocity error
        tmp<volVectorField> Uerr_;

        //- Pressure error
        tmp<volScalarField> pErr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        simpleSolver(const simpleSolver&);

        //- Disallow default bitwise assignment
        void operator=(const simpleSolver&);


public:

    // Constructors

        //- Construct for fvMesh
        simpleSolver(const fvMesh& mesh);


    //- Destructor
    ~simpleSolver();


    // Member Functions

        //- Outer iteration loop control
        bool loop()
        {
            return simple_.loop();
        }

        //- Return the velocity matrix
        tmp<fvVectorMatrix> UEqn(const bool addErr=true);

        //- Solve the system of equations for the specified number of iterations
        void solve(const int nIter);

        //- Calculate and return the velocity residual
        tmp<volVectorField> Ures(const bool addErr=true);

        //- Calculate and return the pressure residual
        tmp<volScalarField> pRes(const bool addErr=true);

        //- Restrict the fine-mesh residuals and set the coarse-mesh errors
        void setError
        (
            simpleSolver& fineEqns,
            const meshToMesh& mapper
        );

        //- Calculate and return the velocity correction
        tmp<volVectorField> Ucorr() const
        {
            return U_ - U0_();
        }

        //- Calculate and return the pressure correction
        tmp<volScalarField> pcorr() const
        {
            return p_ - p0_();
        }

        tmp<volScalarField> kcorr() const
        {
            return mesh_.lookupObject<volScalarField>("k") - k0_();
        }

        tmp<volScalarField> omegacorr() const
        {
            return mesh_.lookupObject<volScalarField>("omega") - omega0_();
        }

        tmp<volScalarField> nutcorr() const
        {
            return mesh_.lookupObject<volScalarField>("nut") - nut0_();
        }

        //- Prolongate the coarse-mesh correction to the fine-mesh
        void correct(const simpleSolver& coarseEqns, const meshToMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
