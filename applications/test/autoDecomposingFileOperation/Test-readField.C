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

#include "argList.H"
#include "volFields.H"
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

    typedef GeometricField<scalar, unallocatedFvPatchField, unallocatedFvMesh>
    uVolScalarField;

    IOobject io
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );


    // Read the undecomposed boundary


    //typedef PtrList<labelList> unallocatedFvBoundaryMesh;
    const unallocatedFvBoundaryMesh boundary;

    unallocatedFvMesh uMesh
    (
        mesh,
        mesh.thisDb(),
        100,
        boundary,
        mesh.globalData()
    );
//    uVolScalarField(io, uMesh);

    volScalarField p
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

    return 0;
}


// ************************************************************************* //
