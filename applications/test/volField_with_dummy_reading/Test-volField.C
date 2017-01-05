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

Application
    Test-volField

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "IOstreams.H"
#include "argList.H"
#include "Time.H"
#include "loadOrCreateMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    //#include "createMesh.H"
    autoPtr<fvMesh> meshPtr
    (
        loadOrCreateMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.constant(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        )
    );
    fvMesh& mesh = meshPtr();


//    {
//        OStringStream os;
//        {
//            autoPtr<fvMesh> dummyMeshPtr
//            (
//                volMesh::New
//                (
//                    IOobject
//                    (
//                        fvMesh::defaultRegion,
//                        runTime.constant(),
//                        runTime,
//                        Foam::IOobject::MUST_READ
//                    ),
//                    mesh,
//                    true
//                )
//            );
//
//            volMesh::write(os, dummyMeshPtr());
//        }
//        {
//            IStringStream is(os.str());
//            autoPtr<fvMesh> dummyMeshPtr
//            (
//                volMesh::New
//                (
//                    IOobject
//                    (
//                        fvMesh::defaultRegion,
//                        runTime.constant(),
//                        runTime,
//                        Foam::IOobject::MUST_READ
//                    ),
//                    is,
//                    true
//                )
//            );
//        }
//    }





    Pout<< "Reading field p\n" << endl;
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

    DebugVar(p);

    return 0;
}


// ************************************************************************* //
