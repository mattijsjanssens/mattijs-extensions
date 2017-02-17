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
    Test-meshToMeshInterpolation

Description
    Use interpolative mapping.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "directMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "map volume fields from one mesh to another"
    );
    argList::noParallel();
    argList::validArgs.append("sourceCase");
    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath = args[1];
    const fileName rootDirSource = casePath.path().toAbsolute();
    const fileName caseDirSource = casePath.name();

    #include "createTimes.H"
    #include "setTimeIndex.H"

    Info<< "Create meshes\n" << endl;

    fvMesh meshSource
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTimeSource.timeName(),
            runTimeSource,
            IOobject::MUST_READ
        )
    );

    fvMesh meshTarget
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTimeTarget.timeName(),
            runTimeTarget,
            IOobject::MUST_READ
        )
    );


    Info<< "Source mesh size: " << meshSource.nCells() << tab
        << "Target mesh size: " << meshTarget.nCells() << nl << endl;


    directMethod meshToMesh(meshSource, meshTarget);

    labelListList srcToTgtAddr;
    scalarListList srcToTgtWght;
    labelListList tgtToSrcAddr;
    scalarListList tgtToSrcWght;
    meshToMesh.calculate
    (
        srcToTgtAddr,
        srcToTgtWght,
        tgtToSrcAddr,
        tgtToSrcWght
    );

    tgtToSrcAddr.clear();
    tgtToSrcWght.clear();

    const pointField& srcCc = meshSource.cellCentres();
    const pointField& tgtCc = meshTarget.cellCentres();


    //DebugVar(srcToTgtAddr);
    //DebugVar(srcToTgtWght);

    volTensorField gradCc(fvc::grad(1.0*meshTarget.C()));

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

            label tgtCelli = tgtCells[0];
            const point& tgtPt = tgtCc[tgtCelli];

            point interpolatedCc = srcPt + ((tgtPt-srcPt)&gradCc[srcCelli]);
            Pout<< "tgt cc:" << interpolatedCc << endl;
        }
        else
        {
            Pout<< "srcCelli:" << srcPt << " does not overlap" << endl;
        }
    }


    //DebugVar(tgtToSrcAddr);
    //DebugVar(tgtToSrcWght);

//     forAll(tgtToSrcAddr, tgtCelli)
//     {
//         const labelList& srcCells = tgtToSrcAddr[tgtCelli];
//
//         if (srcCells.size())
//         {
//             Pout<< "tgtcell:" << tgtCc[tgtCelli]
//                 << " overlaps with src cells:" << pointField(srcCc, srcCells)
//                 << endl;
//         }
//         else
//         {
//             Pout<< "tgtcell:" << tgtCc[tgtCelli]
//                 << " does not overlap" << endl;
//         }
//     }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
