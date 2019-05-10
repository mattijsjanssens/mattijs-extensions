/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    testCompactIOList

Description
    Simple demonstration and test application for the CompactIOList container

\*---------------------------------------------------------------------------*/

#include "ListCompactIO.H"
#include "CompactIOUList.H"
#include "IOstreams.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef CompactIOUList<face, label> faceCompactIOUList;

namespace Foam
{
    defineTemplateTypeNameAndDebugWithName
    (
        faceCompactIOUList,
        "faceCompactIOUList",
        0
    );
}


//  Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
// 
//     IOstream::streamFormat format=IOstream::BINARY;
//     // IOstream::streamFormat format=IOstream::ASCII;

    faceList fcs(2);
    fcs[0] = face(identity(3));
    fcs[1] = face(identity(4));

DebugVar(fcs);

    ListCompactIO<face, label> lst(fcs);

    const face& f = lst[0];
    DebugVar(f);

    {
        // Read
        faceCompactIOUList faces3
        (
            IOobject
            (
                "faces",
                runTime.constant(),
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        forAll(faces3, i)
        {
            DebugVar(faces3[i]);
        }

        Info<< "Read new format " << faces3.size() << " faceList in = "
            << runTime.cpuTimeIncrement() << " s" << nl << endl;
    }

    return 0;
}


// ************************************************************************* //
