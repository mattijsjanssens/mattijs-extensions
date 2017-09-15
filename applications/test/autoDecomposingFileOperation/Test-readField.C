/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Test-readField

Description
    Read volScalarField

\*---------------------------------------------------------------------------*/

#include "GeometricField.H"
#include "argList.H"
#include "uVolFields.H"
#include "unallocatedFvBoundaryMesh.H"
#include "unallocatedFvMesh.H"
#include "unallocatedFvPatchField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOobject io
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );


    // Read the undecomposed boundary
    const fvBoundaryMesh& fp = mesh.boundary();

    // Construct an unallocated boundary. (the references to the
    // original mesh are only used for naming)
    unallocatedFvBoundaryMesh boundary;
    boundary.setSize(3);
    boundary.set(0, new unallocatedGenericFvPatch(fp[0].patch(), 3, fp));
    boundary.set(1, new unallocatedGenericFvPatch(fp[1].patch(), 9, fp));
    boundary.set(2, new unallocatedGenericFvPatch(fp[2].patch(), 18, fp));


    // Construct the mesh
    unallocatedFvMesh uMesh
    (
        mesh,
        mesh.thisDb(),
        9,                      // nCells
        boundary,
        mesh.globalData()
    );


    uVolScalarField p(io, uMesh);
    DebugVar(p);

//    volScalarField p
//    (
//        IOobject
//        (
//            "p",
//            runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );
//    volVectorField U
//    (
//        IOobject
//        (
//            "U",
//            runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );

    return 0;
}


// ************************************************************************* //
