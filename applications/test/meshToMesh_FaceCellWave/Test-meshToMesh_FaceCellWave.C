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
    globalMeshDataTest

Description
    Test global point communication

\*---------------------------------------------------------------------------*/

#include "globalMeshData.H"
#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "mapDistribute.H"
//#include "meshToMeshData.H"
//#include "FaceCellWave.H"
#include "waveMethod.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("sourceCase");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"


    fileName casePath = args[1];
    const fileName rootDirSource = casePath.path().toAbsolute();
    const fileName caseDirSource = casePath.name();
    DebugVar(rootDirSource);
    DebugVar(caseDirSource);
    Time runTimeSource
    (
        Time::controlDictName,
        rootDirSource,
        caseDirSource
    );
    polyMesh meshSource
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTimeSource.timeName(),
            runTimeSource,
            IOobject::MUST_READ
        )
    );


    waveMethod calcMethod(meshSource, mesh);

    labelListList srcToTgtAddr;
    scalarListList srcToTgtWght;
    labelListList tgtToTgtAddr;
    scalarListList tgtToTgtWght;
    calcMethod.calculate
    (
        srcToTgtAddr,
        srcToTgtWght,
        tgtToTgtAddr,
        tgtToTgtWght
    );

    const pointField& srcCc = meshSource.cellCentres();
    const pointField& tgtCc = mesh.cellCentres();

    forAll(srcToTgtAddr, srcCelli)
    {
        const point& srcPt = srcCc[srcCelli];

        const labelList& tgtCells = srcToTgtAddr[srcCelli];

        if (tgtCells.size())
        {
            Pout<< "srcCelli:" << srcPt
                << " is inside tgt cells:"
                << pointField(tgtCc, tgtCells)
                << endl;
        }
        else
        {
            Pout<< "srcCelli:" << srcPt << " does not overlap" << endl;
        }
    }

    forAll(tgtToTgtAddr, tgtCelli)
    {
        const point& tgtPt = tgtCc[tgtCelli];

        const labelList& srcCells = tgtToTgtAddr[tgtCelli];

        if (srcCells.size())
        {
            Pout<< "tgtCelli:" << tgtPt
                << " is inside tgt cells:"
                << pointField(srcCc, srcCells)
                << endl;
        }
        else
        {
            Pout<< "tgtCelli:" << tgtPt << " does not overlap" << endl;
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
