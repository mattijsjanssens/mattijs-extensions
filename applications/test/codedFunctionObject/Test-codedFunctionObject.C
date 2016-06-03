/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "argList.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "stringOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

//     wordHashSet vars;
//     {
//         string a("// SHA1 = ${SHA1sum}");
//         stringOps::variables(vars, a);
//         DebugVar(vars);
//     }
//     {
//         string b("const char* const ${typeName}FixedValuePointPatch${FieldType}::SHA1sum =");
//         vars.clear();
//         stringOps::variables(vars, b);
//         DebugVar(vars);
//     }

    const pointMesh& pMesh = pointMesh::New(mesh);

    pointScalarField pointField
    (
        IOobject
        (
            "pointField",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        pMesh
    );

DebugVar(pointField);


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
