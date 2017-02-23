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

#include "meshToMesh.H"
#include "fvCFD.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//autoPtr<mapDistribute>
void calcGlobalCellCells
(
    const globalIndex& globalNumbering,
    const fvMesh& mesh,
    const PackedBoolList& needStencil,
    labelListList& cellCells
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const cellList& cells = mesh.cells();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList globalCellIDs(mesh.nCells());
    forAll(globalCellIDs, celli)
    {
        globalCellIDs[celli] = globalNumbering.toGlobal(celli);
    }

    labelList nbrGlobalCellIDs;
    syncTools::swapBoundaryCellList
    (
        mesh,
        globalCellIDs,
        nbrGlobalCellIDs
    );

    cellCells.setSize(mesh.nCells());
    forAll(cellCells, celli)
    {
        if (needStencil[celli])
        {
            const cell& cFaces = cells[celli];

            labelList& stencil = cellCells[celli];
            stencil.setSize(1+cFaces.size());
            label index = 0;
            stencil[index++] = globalNumbering.toGlobal(celli);
            forAll(cFaces, i)
            {
                label facei = cFaces[i];
                if (mesh.isInternalFace(facei))
                {
                    label nbrCelli = own[facei];
                    if (nbrCelli == celli)
                    {
                        nbrCelli = nei[facei];
                    }
                    stencil[index++] = globalNumbering.toGlobal(nbrCelli);
                }
                else
                {
                    label bFacei = facei-mesh.nInternalFaces();
                    label patchi = pbm.patchID()[bFacei];
                    if (pbm[patchi].coupled())
                    {
                        stencil[index++] = nbrGlobalCellIDs[bFacei];
                    }
                }
            }
            stencil.setSize(index);
        }
    }

//    List<Map<label>> compactMap;
//    return autoPtr<mapDistribute>
//    (
//        new mapDistribute
//        (
//            globalNumbering,
//            cellCells,
//            compactMap
//        )
//    );
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "map volume fields from one mesh to another"
    );
    argList::validArgs.append("sourceCase");
    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.caseName());

    fileName casePath = args[1];
    const fileName rootDirSource = casePath.path().toAbsolute();
    const word myDir(word("processor") + Foam::name(Pstream::myProcNo()));
    const fileName caseDirSource = casePath.name()/myDir;

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


    Pout<< "Source mesh size: " << meshSource.nCells() << tab
        << "Target mesh size: " << meshTarget.nCells() << nl << endl;

    // Ignore all patches by making them cuttingPatches.
    wordList cuttingPatches(meshTarget.boundaryMesh().names());

    meshToMesh mapper
    (
        meshSource,
        meshTarget,
        meshToMesh::imDirect,
        HashTable<word>(),
        cuttingPatches
    );

    const labelListList& srcToTgtAddr = mapper.srcToTgtCellAddr();
    const pointField& srcCc = meshSource.cellCentres();
    const pointField& tgtCc = meshTarget.cellCentres();


    //DebugVar(srcToTgtAddr);
    //DebugVar(srcToTgtWght);

    volTensorField gradCc(fvc::grad(1.0*meshTarget.C()));

DebugVar(mapper.singleMeshProc());

    if (mapper.singleMeshProc() != -1)
    {
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
    }
    else
    {
        // Get target cells local
        pointField tgtCc(meshTarget.cellCentres());
        mapper.tgtMap().distribute(tgtCc);


        pointField samples(srcCc.size(), vector::max);

        forAll(srcToTgtAddr, srcCelli)
        {
            const point& srcPt = srcCc[srcCelli];

            const labelList& slots = srcToTgtAddr[srcCelli];

            if (slots.size())
            {
                Pout<< "srcCelli:" << srcPt
                    << " is inside tgt cells:"
                    << pointField(tgtCc, slots)
                    << endl;

                label sloti = slots[0];
                const point& tgtPt = tgtCc[sloti];

                point interpolatedCc = srcPt + ((tgtPt-srcPt)&gradCc[srcCelli]);
                Pout<< "tgt cc:" << interpolatedCc << endl;

                samples[sloti] = srcPt;
            }
            else
            {
                Pout<< "srcCelli:" << srcPt << " does not overlap" << endl;
            }
        }

        mapper.tgtMap().reverseDistribute(meshTarget.nCells(), samples);

        PackedBoolList needStencil(meshTarget.nCells());
        forAll(samples, tgtCelli)
        {
            if (samples[tgtCelli] != vector::max)
            {
                Pout<< "** target cell:" << meshTarget.cellCentres()[tgtCelli]
                    << " contains source cellcc:" << samples[tgtCelli]
                    << endl;
                needStencil[tgtCelli] = true;
            }
        }

        // Determine stencil (in global indices) for selected cells

        globalIndex globalTgtNumbers(meshTarget.nCells());

        //XXXXX
        labelListList tgtCellCells;
        calcGlobalCellCells
        (
            globalTgtNumbers,
            meshTarget,
            needStencil,
            tgtCellCells
        );

        // Pull back to src
        mapper.tgtMap().distribute(tgtCellCells);

        // Src determine cell centres of stencil and weighting factors
        labelListList srcCellToTgtCells(meshSource.nCells());


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
