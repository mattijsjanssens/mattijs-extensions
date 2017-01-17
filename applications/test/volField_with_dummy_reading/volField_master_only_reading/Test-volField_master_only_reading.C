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
    Test-volField

Description
    The idea is that argList allocates a fileServer (through command-line
    arguments or environment vars). This has operations to do
    - file existence checking
    - mkDir, rm etc.
    - open an IFstream
    - open an OFstream
    and everywhere instead of constructing an I/OFstream we do a
        fileServer::NewIFstream(io)
        fileServer::NewOFstream(io)

\*---------------------------------------------------------------------------*/

//#include "IOstreams.H"
#include "argList.H"
#include "Time.H"
//#include "volFields.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "OFstream.H"
//#include "masterOFstream.H"
#include "fileOperation.H"
#include "localFileOperation.H"
#include "masterFileOperation.H"
#include "pointFields.H"

using namespace Foam;

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOList<List<char>>, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

Pout<< "std::streamoff:" << sizeof(std::streamoff) << endl;

    // Test IOList writing
    {
        IOList<List<char>> bufs
        (
            IOobject
            (
                "bufs",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            2
        );

        forAll(bufs, proci)
        {
            List<char>& buf = bufs[proci];

            buf.setSize(10);
            forAll(buf, i)
            {
                buf[i] = proci+'a'+i;
            }
        }

        bufs.write();

        return 0;
    }



    const pointMesh& pm = pointMesh::New(mesh);

    IOobject io
    (
        "pointDisplacement",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    pointVectorField pointDisplacement(io, pm);

DebugVar(pointDisplacement);

    IOdictionary sol
    (
        IOobject
        (
            "fvSolution",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
DebugVar(sol);

Pout<< "*** WRITING **" << endl;

        runTime++;
        pointDisplacement.instance() = runTime.timeName();
        pointDisplacement.regIOobject::write();

Pout<< "*** READING **" << endl;

    {
        pointVectorField pointDisplacement2
        (
            IOobject
            (
                "pointDisplacement",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            pm
        );
        DebugVar(pointDisplacement2);
    }

    return 0;
}


// ************************************************************************* //
