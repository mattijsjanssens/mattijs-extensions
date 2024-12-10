/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 M. Janssens
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

#include "decompositionMethod.H"
#include "decompositionGAMGAgglomeration.H"
//#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "decompositionMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        decompositionGAMGAgglomeration,
        lduMesh
    );

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        decompositionGAMGAgglomeration,
        geometry
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionGAMGAgglomeration::decompositionGAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    GAMGAgglomeration(mesh, controlDict),
    decomposerPtr_
    (
        decompositionMethod::New
        (
            controlDict.optionalSubDict(type() + "Coeffs")
        )
    )   // regionName
{
    DebugVar(controlDict);

    agglomerate
    (
        nCellsInCoarsestLevel_,
        0,          //mesh,
        scalarField(mesh.lduAddr().upperAddr().size(), scalar(1.0)),
        true
    );
}


Foam::decompositionGAMGAgglomeration::decompositionGAMGAgglomeration
(
    const lduMesh& mesh,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
:
    decompositionGAMGAgglomeration(mesh, controlDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::decompositionGAMGAgglomeration::localCellCells
(
    const lduAddressing& addr,
    const bitSet& isBlockedFace,
    CompactListList<label>& cellCells
)
{
    // Since have only internal addressing cannot get global cell connectivity?

    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

    labelList nNbrs(addr.size(), Zero);
    forAll(nbr, facei)
    {
        if (!isBlockedFace[facei])
        {
            nNbrs[nbr[facei]]++;
            nNbrs[own[facei]]++;
        }
    }
    // Calculate&store offsets
    cellCells.resize_nocopy(nNbrs);

    auto& values = cellCells.values();
    auto& offsets = cellCells.offsets();

    // Fill in neighbours
    nNbrs = Zero;
    forAll(nbr, facei)
    {
        if (!isBlockedFace[facei])
        {
            const label n = nbr[facei];
            const label o = own[facei];
            values[offsets[n]+(nNbrs[n]++)] = o;
            values[offsets[o]+(nNbrs[o]++)] = n;
        }
    }
}
void Foam::decompositionGAMGAgglomeration::localCellCells
(
    const lduAddressing& addr,
    const labelList& regions,   // marker
    const label regioni,        // which marker to keep
    CompactListList<label>& cellCells
)
{
    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

    labelList oldToNew(regions.size(), -1);
    label nCells = 0;
    forAll(regions, i)
    {
        if (regions[i] == regioni)
        {
            oldToNew[i] = nCells++;
        }
    }

    labelList nNbrs(nCells, Zero);
    forAll(nbr, facei)
    {
        const label n = oldToNew[nbr[facei]];
        const label o = oldToNew[own[facei]];
        if (n != -1 && o != -1)
        {
            nNbrs[n]++;
            nNbrs[o]++;
        }
    }
    // Calculate&store offsets
    cellCells.resize_nocopy(nNbrs);

    auto& values = cellCells.values();
    auto& offsets = cellCells.offsets();

    // Fill in neighbours
    nNbrs = Zero;
    forAll(nbr, facei)
    {
        const label n = oldToNew[nbr[facei]];
        const label o = oldToNew[own[facei]];
        if (n != -1 && o != -1)
        {
            values[offsets[n]+(nNbrs[n]++)] = o;
            values[offsets[o]+(nNbrs[o]++)] = n;
        }
    }
}


Foam::tmp<Foam::labelField> Foam::decompositionGAMGAgglomeration::agglomerate
(
    label& nCoarseCells,
    const lduAddressing& addr,
    const scalarField& faceWeights
) const
{
    const labelUList& nbr = addr.upperAddr();
    //const labelUList& own = addr.lowerAddr();

    const globalIndex globalCells
    (
        addr.size(),
        UPstream::worldComm,
        false           // local only
    );


    nCoarseCells = addr.size()/2;
Pout<< "Agglomerating from " << addr.size() << " to " << nCoarseCells
    << endl;

    // Construct globalCellCells for calling decomposition method.
    CompactListList<label> globalCellCells;
    localCellCells(addr, bitSet(nbr.size()), globalCellCells);

    //for
    //(
    //    label celli = 0;
    //    celli < globalCellCells.offsets().size()-1;
    //    celli++
    //)
    //{
    //    Pout<< "For cell:" << celli
    //        << " nbrs:" << flatOutput(globalCellCells[celli])
    //        << endl;
    //}


    tmp<labelField> tfld(new labelField());

    decompositionMethod& decomposer = decomposerPtr_();
    decomposer.nDomains(nCoarseCells);
    tfld.ref() = decomposer.decompose
    (
        globalCellCells,
        pointField(globalCellCells.size(), Zero)
    );

    return tfld;
}

//void Foam::decompositionGAMGAgglomeration::agglomerate
//(
//    const lduAddressing& startMesh,
//    const scalarField& startFaceWeights,
//
//    PtrList<labelList>& agglomeration
//)
//{
//    decompositionMethod& decomposer = decomposerPtr_();
//
//    agglomeration.resize_nocopy(maxLevels_);
//
//
//    const globalIndex globalCells
//    (
//        addr.size(),
//        UPstream::worldComm,
//        false           // local only
//    );
//
//
//    // Start geometric agglomeration from the given faceWeights
//    scalarField faceWeights = startFaceWeights;
//
//    CompactListList<label> cellCells;
//    localCellCells(addr, cellCells);
//
//    // Decompose into coarsest level
//    decomposer.nDomains(nCoarseCells);
//    agglomeration.set
//    (
//        0,
//        decomposer.decompose
//        (
//            cellCells,
//            pointField(cellCells.size(), Zero)
//        )
//    );
//    const auto& agglom = agglomeration[0];
//
//    //// Do each resulting domain in turn (inefficient)
//    //const label nDecomp = max(agglom);
//    //for (label i = 0; i < nDecomp; i++)
//    //{
//    //    CompactListList<label> subCellCells;
//    //    localCellCells
//    //    (
//    //        startMesh,
//    //        agglom,
//    //        i,
//    //        subCellCells
//    //    );
//    //
//    //    tmp<labelField> tsubDecomp
//    //    (
//    //        decomposer.decompose
//    //        (
//    //            subCellCells,
//    //            pointField(subCellCells.size(), Zero)
//    //        )
//    //    );
//    //    const auto& subDecomp = tsubDecomp();
//    //}
//
//
//    // Do decomposition by generating cell-cells without connections
//    const labelUList& nbr = addr.upperAddr();
//    const labelUList& own = addr.lowerAddr();
//    bitSet isBlockedFace(nbr.size());
//    forAll(nbr, facei)
//    {
//        const label n = nbr[facei];
//        const label o = own[facei];
//        if (agglom[n] != agglom[o])
//        {
//            isBlockedFace.set(facei);
//        }
//    }
//
//    CompactListList<label> cellCells;
//    localCellCells
//    (
//        addr,
//        isBlockedFace,
//        cellCells
//    );
//
//    decomposer.nDomains(2);
//    tmp<labelField> tdecomp = decomposer.decompose
//    (
//        cellCells,
//        pointField(cellCells.size(), Zero)
//    );
//
//    forAll(agglom, celli)
//    {
//        agglom[celli] = decomposer.nDomains()*agglom[celli]+tdecomp()[celli];
//    }
//
//    while()
//    {
//
//
//
//    }


void Foam::decompositionGAMGAgglomeration::agglomerate
(
    const label nCellsInCoarsestLevel,
    const label startLevel,
    const scalarField& startFaceWeights,
    const bool doProcessorAgglomerate
)
{
    // Straight copy of pairGAMGAgglomeration::agglomerate without the
    // while loop.

    if (nCells_.size() < maxLevels_)
    {
        // See compactLevels. Make space if not enough
        nCells_.resize(maxLevels_);
        restrictAddressing_.resize(maxLevels_);
        nFaces_.resize(maxLevels_);
        faceRestrictAddressing_.resize(maxLevels_);
        faceFlipMap_.resize(maxLevels_);
        nPatchFaces_.resize(maxLevels_);
        patchFaceRestrictAddressing_.resize(maxLevels_);
        meshLevels_.resize(maxLevels_);
        // Have procCommunicator_ always, even if not procAgglomerating.
        // Use value -1 to indicate nothing is proc-agglomerated
        procCommunicator_.resize(maxLevels_ + 1, -1);
        if (processorAgglomerate())
        {
            procAgglomMap_.resize(maxLevels_);
            agglomProcIDs_.resize(maxLevels_);
            procCommunicator_.resize(maxLevels_);
            procCellOffsets_.resize(maxLevels_);
            procFaceMap_.resize(maxLevels_);
            procBoundaryMap_.resize(maxLevels_);
            procBoundaryFaceMap_.resize(maxLevels_);
        }
    }


    // Start geometric agglomeration from the given faceWeights
    scalarField faceWeights = startFaceWeights;

    // Agglomerate until the required number of cells in the coarsest level
    // is reached
    label nCreatedLevels = startLevel;

    while (nCreatedLevels < maxLevels_ - 1)
    {
        if (!hasMeshLevel(nCreatedLevels))
        {
            FatalErrorInFunction<< "No mesh at nCreatedLevels:"
                << nCreatedLevels
                << exit(FatalError);
        }

        const auto& fineMesh = meshLevel(nCreatedLevels);


        label nCoarseCells = -1;

        // Agglomerate single level : determine per fineMesh cell
        // which cluster it goes to. Sets number of clusters.
        tmp<labelField> finalAgglomPtr = agglomerate
        (
            nCoarseCells,
            fineMesh.lduAddr(),
            faceWeights
        );

        DebugVar(nCoarseCells);
        Pout<< "resulting agglom:" << finalAgglomPtr() << endl;

        if
        (
            continueAgglomerating
            (
                nCellsInCoarsestLevel,
                finalAgglomPtr().size(),
                nCoarseCells,
                fineMesh.comm()
            )
        )
        {
            nCells_[nCreatedLevels] = nCoarseCells;
            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
        }
        else
        {
            break;
        }

        // Create coarse mesh
        agglomerateLduAddressing(nCreatedLevels);

        // Agglomerate the faceWeights field for the next level
        {
            scalarField aggFaceWeights
            (
                meshLevels_[nCreatedLevels].upperAddr().size(),
                0.0
            );

            restrictFaceField
            (
                aggFaceWeights,
                faceWeights,
                nCreatedLevels
            );

            faceWeights = std::move(aggFaceWeights);
        }

        nCreatedLevels++;
    }

    // Shrink the storage of the levels to those created
    compactLevels(nCreatedLevels, doProcessorAgglomerate);
}


// ************************************************************************* //
