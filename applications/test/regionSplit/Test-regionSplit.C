/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Mattijs Janssens
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

#include "fvMesh.H"
#include "volFields.H"
#include "argList.H"
#include "Time.H"
#include "zeroGradientFvPatchFields.H"
#include "regionSplit.H"
#include "faceSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("faceSet");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word setName = args[1];
    faceSet setFaces(mesh, setName);

    boolList isBlocked(mesh.nFaces(), false);
    forAllConstIter(faceSet, setFaces, iter)
    {
        isBlocked[iter.key()] = true;
    }


    regionSplit regions(mesh, isBlocked);

    volScalarField fld
    (
        IOobject
        (
            "cellRegion",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    scalarField& f = fld.ref();
    forAll(regions, celli)
    {
        f[celli] = regions[celli];
    }
    fld.correctBoundaryConditions();
    fld.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
