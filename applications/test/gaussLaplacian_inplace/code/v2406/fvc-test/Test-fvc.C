/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fvcSurfaceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    const volScalarField one
    (
        IOobject
        (
            "one",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<scalar>(dimLength, 1.0),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    DebugVar(one);

    scalarField a(3, 123);
    const scalarField b({1, 2, 3});
    const scalarField c({2, 4, 5});

    multiplySubtract(a, b, c);
DebugVar(a);

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
DebugVar(T);
    const volScalarField T2("T2", T);
DebugVar(T2);


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    const volVectorField U2("U2", U);

    #include "createPhi.H"

    // Grad
    Pout<< "tgradT:" << fvc::grad(T) << endl;
    Pout<< "tgradT2:" << fvc::grad(T2) << endl;

    // Div
    Pout<< "tdivU:" << fvc::div(U) << endl;
    Pout<< "tdivU2:" << fvc::div(U2) << endl;
    // Convection
    Pout<< "tdiv(phi,U):" << fvc::div(phi,U) << endl;
    Pout<< "tdiv(phi,U2):" << fvc::div(phi,U2) << endl;
    // Laplacian
    Pout<< "tlaplacianT:" << fvc::laplacian(T) << endl;
    Pout<< "tlaplacianT2:" << fvc::laplacian(T2) << endl;
    // LeastSquares

return 0;

cpuTime timer;
{
    Pout<< "tgradT:" << endl;   //<< tgradT() << endl;
    for (label i = 0; i < 100; i++)
    {
        tmp<volVectorField> tgradT(fvc::grad(T));
    }
    Pout<< "tgradT:" << timer.cpuTimeIncrement() << endl;
}
{
    Pout<< "tgradT2:" << endl;
    for (label i = 0; i < 100; i++)
    {
        tmp<volVectorField> tgradT2(fvc::grad(T2));
    }
    Pout<< "tgradT2:" << timer.cpuTimeIncrement() << endl;
}

{
    Pout<< "tlaplacianT:" << endl;   //<< tgradT() << endl;
    for (label i = 0; i < 100; i++)
    {
        tmp<volScalarField> tlaplacianT(fvc::laplacian(one, T));
    }
    Pout<< "tlaplacianT:" << timer.cpuTimeIncrement() << endl;
}
{
    Pout<< "tlaplacianT2:" << endl;   //<< tgradT() << endl;
    for (label i = 0; i < 100; i++)
    {
        tmp<volScalarField> tlaplacianT2(fvc::laplacian(one, T2));
    }
    Pout<< "tlaplacianT2:" << timer.cpuTimeIncrement() << endl;
}

{
    Pout<< "tsnGradT:" << endl;   //<< tgradT() << endl;
    for (label i = 0; i < 100; i++)
    {
        tmp<surfaceScalarField> tsnGradT(fvc::snGrad(T));
    }
    Pout<< "tsnGradT:" << timer.cpuTimeIncrement() << endl;
}
{
    Pout<< "tsnGradT2:" << endl;   //<< tgradT() << endl;
    for (label i = 0; i < 100; i++)
    {
        tmp<surfaceScalarField> tsnGradT(fvc::snGrad(T2));
    }
    Pout<< "tsnGradT2:" << timer.cpuTimeIncrement() << endl;
}

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
