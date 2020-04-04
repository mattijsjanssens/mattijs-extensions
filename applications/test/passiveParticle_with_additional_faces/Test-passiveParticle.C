/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    testPassiveParticle

Description
    Test cloud of passive particles.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "passiveParticleCloud.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noFunctionObjects();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word cloudName = "myFirstCloud";

    // Start with empty cloud
    passiveParticleCloud particles
    (
        mesh,
        cloudName,
        IDLList<passiveParticle>()
    );
    Pout<< "Starting particles:" << particles.size() << endl;

    Pout<< "Adding a particle." << endl;
    particles.addParticle(new passiveParticle(mesh, point(0.1, 0.01, 0.001)));

    for (const passiveParticle& p : particles)
    {
        Pout<< "    " << p.position() << " cell:" << p.cell()
            << " origProc:" << p.origProc()
            << " origId:" << p.origId()
            << " coord:" << p.coordinates()
            << " transform:" << p.currentTetTransform()
            << endl;
    }

    runTime.printExecutionTime(Info);

    //++runTime;
    //Pout<< "Writing particles to time " << runTime.timeName() << endl;
    //particles.write();

    const pointField immPoints
    ({
        point(0, 0, 0.99),
        point(1, 0, 0.99),
        point(0, 1, 0.99),
        point(1, 1, 0.99)
    });
    const triFaceList immFaces
    ({
        triFace(0, 1, 2),
        triFace(2, 1, 3)
    });
    labelListList cellTriangles(2);
    cellTriangles[0] = labelList({0, 1});

    {
        OBJstream os(runTime.path()/"cellTriangles.obj");
        Pout<< "Writing immersed triangles to " << os.name() << endl;
        for (const labelList& cTris : cellTriangles)
        {
            for (const label triI : cTris)
            {
                os.write(immFaces[triI].tri(immPoints), false);
            }
        }
    }


    // Track the particles to the next face
    for (passiveParticle& p : particles)
    {
        const point start(p.position());
        const point end(0.1, 0.01, 1.03);

        Pout<< "Tracking " << start << " cell:" << p.cell()
            << " origProc:" << p.origProc()
            << " origId:" << p.origId()
            << " to:" << end
            << " coord:" << p.coordinates()
            << " transform:" << p.currentTetTransform()
            << endl;

        const scalar s = p.trackToFace
        (
            end-start,
            1.0,
            immPoints,
            immFaces,
            cellTriangles[p.cell()]
        );

        Pout<< "    After:" << " cell:" << p.cell()
            << " pos:" << p.position()
            << " fraction:" << s
            << " coord:" << p.coordinates()
            << " transform:" << p.currentTetTransform()
            << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
