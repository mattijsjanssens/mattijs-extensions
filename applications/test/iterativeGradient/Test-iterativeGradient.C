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

Application
    sphericalTensorFieldTest

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvCFD.H"
#include "skewCorrectionVectors.H"
//#include "volFields.H"
//#include "surfaceFields.H"
#include "pisoControl.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field p\n" << endl;
    volScalarField vfld
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    vfld.dimensions().reset(dimLength);
    vfld == mag(mesh.C());


{
    volVectorField gradP
    (
        IOobject
        (
            "gradP",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(vfld)
    );
}



//DebugVar(vfld);
//    const skewCorrectionVectors& skv = skewCorrectionVectors::New(mesh);
//
////DebugVar(skv());
//
//    surfaceScalarField sfld
//    (
//        IOobject
//        (
//            "sfld",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        linearInterpolate(vfld)
//    );
//    volVectorField gradP
//    (
//        IOobject
//        (
//            "gradP",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        fvc::grad(sfld)
//    );
//
//    gradP.write();
//
//    pisoControl piso(mesh);
//
//    Info<< "\nStarting time loop\n" << endl;
//
//    while (runTime.loop())
//    {
//        Info<< "Time = " << runTime.timeName() << nl << endl;
//        surfaceVectorField sgradP("sgradP", linearInterpolate(gradP));
//DebugVar(linearInterpolate(vfld));
//        surfaceScalarField corr("skewCorr", (skv()&sgradP));
//        corr.dimensions().reset(vfld.dimensions());
//        sfld = linearInterpolate(vfld)+corr;
//DebugVar(sfld);
//        gradP = fvc::grad(sfld);
//
//        runTime.write();
//DebugVar(gradP);
//    }

    return 0;
}


// ************************************************************************* //
