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

Foam::label Foam::decompositionGAMGAgglomeration::regionSplit
(
    const lduAddressing& addr,
    const bitSet& isBlockedFace,
    labelList& cellRegion
) const
{
    // Determine cell-to-cell addressing, avoiding crossing blocked faces)
    CompactListList<label> cellCells;
    localCellCells(addr, isBlockedFace, cellCells);

    // Walk, assign regions. Taken from regionSplit
    cellRegion.setSize(addr.size(), -1);

    // Start with region 0
    label nLocalRegions = 0;

    for (label seedCelli = 0; seedCelli < cellRegion.size(); ++seedCelli)
    {
        // Find next unset cell - use as seed

        for (; seedCelli < cellRegion.size(); ++seedCelli)
        {
            if (cellRegion[seedCelli] == -1)
            {
                break;
            }
        }

        if (seedCelli >= cellRegion.size())
        {
            break;
        }

        // Seed cell
        cellRegion[seedCelli] = nLocalRegions;


        // Walk nLocalRegions out
        DynamicList<label> changedCells;
        changedCells.append(seedCelli);

        while (changedCells.size())
        {
            DynamicList<label> newChangedCells(changedCells.size());
            for (const label celli : changedCells)
            {
                for (const label nbrCelli : cellCells[celli])
                {
                    if (cellRegion[nbrCelli] == -1)
                    {
                        cellRegion[nbrCelli] = nLocalRegions;
                        newChangedCells.append(nbrCelli);
                    }
                }
            }

            changedCells = newChangedCells;
        }

        ++nLocalRegions; // Next region
    }

    return nLocalRegions;
}
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

    // Assign cellToRegion


}
Foam::labelList Foam::decompositionGAMGAgglomeration::localCellCells
(
    const lduAddressing& addr,
    const labelList& regions,   // marker per cell
    const label regioni,        // which marker to keep
    CompactListList<label>& cellCells
)
{
    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

    if (regions.size() != addr.size())
    {
        FatalErrorInFunction<< "Wrong size" << exit(FatalError);
    }

    // Compact cell numbering for region cells
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
    return oldToNew;
}
Foam::bitSet Foam::decompositionGAMGAgglomeration::blockedFaces
(
    const lduAddressing& addr,
    const labelUList& region        // region per cell
)
{
    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

    bitSet isBlockedFace(nbr.size());
    forAll(nbr, facei)
    {
        isBlockedFace[facei] = (region[own[facei]] != region[nbr[facei]]);
    }
    return isBlockedFace;
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


    nCoarseCells = addr.size()/4;
Pout<< "Agglomerating from " << addr.size() << " to " << nCoarseCells
    << endl;

    // Construct globalCellCells for calling decomposition method.
    CompactListList<label> globalCellCells;
    localCellCells(addr, bitSet(nbr.size()), globalCellCells);

    for
    (
        label celli = 0;
        celli < globalCellCells.offsets().size()-1;
        celli++
    )
    {
        Pout<< "For cell:" << celli
            << " nbrs:" << flatOutput(globalCellCells[celli])
            << endl;
    }


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


//void Foam::decompositionGAMGAgglomeration::agglomerate
//(
//    const label nCellsInCoarsestLevel,
//    const label startLevel,
//    const scalarField& startFaceWeights,
//    const bool doProcessorAgglomerate
//)
//{
//    // Straight copy of pairGAMGAgglomeration::agglomerate without the
//    // while loop.
//
//    if (nCells_.size() < maxLevels_)
//    {
//        // See compactLevels. Make space if not enough
//        nCells_.resize(maxLevels_);
//        restrictAddressing_.resize(maxLevels_);
//        nFaces_.resize(maxLevels_);
//        faceRestrictAddressing_.resize(maxLevels_);
//        faceFlipMap_.resize(maxLevels_);
//        nPatchFaces_.resize(maxLevels_);
//        patchFaceRestrictAddressing_.resize(maxLevels_);
//        meshLevels_.resize(maxLevels_);
//        // Have procCommunicator_ always, even if not procAgglomerating.
//        // Use value -1 to indicate nothing is proc-agglomerated
//        procCommunicator_.resize(maxLevels_ + 1, -1);
//        if (processorAgglomerate())
//        {
//            procAgglomMap_.resize(maxLevels_);
//            agglomProcIDs_.resize(maxLevels_);
//            procCommunicator_.resize(maxLevels_);
//            procCellOffsets_.resize(maxLevels_);
//            procFaceMap_.resize(maxLevels_);
//            procBoundaryMap_.resize(maxLevels_);
//            procBoundaryFaceMap_.resize(maxLevels_);
//        }
//    }
//
//    decompositionMethod& decomposer = decomposerPtr_();
//
//    // Start geometric agglomeration from the given faceWeights
//    scalarField faceWeights = startFaceWeights;
//
//    // Agglomerate until the required number of cells in the coarsest level
//    // is reached
//    label nCreatedLevels = startLevel;
//
//    // Per level the agglomeration
//    DynamicList<labelList> agglomeration(maxLevels_);
//    // Per level the size (= max of agglomeration)
//    DynamicList<label> nCells(maxLevels_);
//
//    const auto& startMesh = meshLevel(startLevel);
//    const lduAddressing& startAddr = startMesh.lduAddr();
//
//    // Note: per level (starting from coarsest), per start mesh cell the
//    //       agglomeration (= region)
//    DynamicList<labelList> startAgglomeration(maxLevels_);
//
//    bitSet isBlockedFace(startAddr.lowerAddr().size());
//
//
//    // Agglomerate current level
//    {
//        // All cells in same agglomeration
//        labelList oldAgglom(startAddr.size(), 0);
//
//        // Calculate cellCells (or rather region-regions)
//        CompactListList<label> cellCells;
//        localCellCells
//        (
//            startAddr,
//            isBlockedFace, // dummy since in same region
//            cellCells
//        );
//        //labelList cellRegion;
//        //regionSplit(startAddr, isBlockedFace, cellRegion);
//
//        decomposer.nDomains(nCellsInCoarsestLevel);
//        const bool oldParRun = UPstream::parRun(false);
//        const labelList cellToProc
//        (
//            decomposer.decompose
//            (
//                cellCells,
//                pointField(cellCells.size(), Zero)
//            )
//        );
//        UPstream::parRun(oldParRun);
//
//        agglomeration.append(cellToProc);
//
//        agglomeration.append(labelList(startAddr.size()));
//        auto& agglom = agglomeration.last();
//        UIndirectList<label>(agglom, oldAgglom) = cellToProc;
//
//        nCells_[nCreatedLevels] = max(cellToProc);
//        restrictAddressing_.set(nCreatedLevels, new labelField(agglom));
//
//        nCreatedLevels++;
//    }
//
///*
//    while (nCreatedLevels < maxLevels_ - 1)
//    {
//        labelList& oldAgglom = agglomeration.last();
//
//        // Calculate cellCells, as if split according to oldAgglom
//        CompactListList<label> cellCells;
//        localCellCells
//        (
//            startAddr,
//            blockedFaces(startAddr, oldAgglom),
//            cellCells
//        );
//
//        if (cellCells.size() > fineMesh.nCells()/4)
//        {
//            break;
//        }
//
//        nCells.append(cellCells.size());
//
//        decomposer.nDomains(4);
//        const bool oldParRun = UPstream::parRun(false);
//        const labelList cellToProc
//        (
//            decomposer.decompose
//            (
//                cellCells,
//                pointField(cellCells.size(), Zero)
//            )
//        );
//        UPstream::parRun(oldParRun);
//
//
//        agglomeration.append(labelList::New(startMesh.size()));
//        auto& agglom = agglomeration.last();
//        UIndirectList<label>(agglom, oldAgglom) = cellToProc;
//
//        nCreatedLevels++;
//    }
//
//
//    // Use agglomeration to set up addressing
//    
//
//
//
//    while (nCreatedLevels < maxLevels_ - 1)
//    {
//        if (!hasMeshLevel(nCreatedLevels))
//        {
//            FatalErrorInFunction<< "No mesh at nCreatedLevels:"
//                << nCreatedLevels
//                << exit(FatalError);
//        }
//
//        const auto& fineMesh = meshLevel(nCreatedLevels);
//
//        {
//            Pout<< "At level:" << nCreatedLevels << endl;
//
//            const auto& addr = fineMesh.lduAddr();
//            const labelUList& nbr = addr.upperAddr();
//            const labelUList& own = addr.lowerAddr();
//            Pout<< "nCells:" << addr.size() << endl;
//            Pout<< "nFaces:" << nbr.size() << endl;
//            forAll(nbr, facei)
//            {
//                Pout<< "    face:" << facei
//                    << " between:" << own[facei]
//                    << " and:" << nbr[facei]
//                    << endl;
//            }
//        }
//
//
//
//
//        label nCoarseCells = -1;
//
//        // Agglomerate single level : determine per fineMesh cell
//        // which cluster it goes to. Sets number of clusters.
//        tmp<labelField> finalAgglomPtr = agglomerate
//        (
//            nCoarseCells,
//            fineMesh.lduAddr(),
//            faceWeights
//        );
//
//        DebugVar(nCoarseCells);
//        Pout<< "resulting agglom:" << finalAgglomPtr() << endl;
//
//        if
//        (
//            continueAgglomerating
//            (
//                nCellsInCoarsestLevel,
//                finalAgglomPtr().size(),
//                nCoarseCells,
//                fineMesh.comm()
//            )
//        )
//        {
//            nCells_[nCreatedLevels] = nCoarseCells;
//            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
//        }
//        else
//        {
//            break;
//        }
//
//        // Create coarse mesh
//        agglomerateLduAddressing(nCreatedLevels);
//
//        // Agglomerate the faceWeights field for the next level
//        {
//            scalarField aggFaceWeights
//            (
//                meshLevels_[nCreatedLevels].upperAddr().size(),
//                0.0
//            );
//
//            restrictFaceField
//            (
//                aggFaceWeights,
//                faceWeights,
//                nCreatedLevels
//            );
//
//            faceWeights = std::move(aggFaceWeights);
//        }
//
//        nCreatedLevels++;
//    }
//
//*/
//
//    // Shrink the storage of the levels to those created
//    compactLevels(nCreatedLevels, doProcessorAgglomerate);
//}

//XXXXX
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

    decompositionMethod& decomposer = decomposerPtr_();

    // Agglomerate until the required number of cells in the coarsest level
    // is reached
    label nCreatedLevels = startLevel;

    // Per level the agglomeration
    DynamicList<labelList> agglomeration(maxLevels_);
    // Per level the size (= max of agglomeration)
    DynamicList<label> nCells(maxLevels_);

    // Current mesh
    autoPtr<lduPrimitiveMesh> levelMeshPtr;

    while (nCreatedLevels < maxLevels_ - 1)
    {
        const lduMesh& mesh =
        (
            nCreatedLevels == startLevel
          ? meshLevel(nCreatedLevels)
          : levelMeshPtr()
        );

        if (mesh.lduAddr().size() < nCellsInCoarsestLevel)
        {
            break;
        }

        labelList dummyAgglom;
        if (nCreatedLevels == startLevel)
        {
            dummyAgglom.setSize(mesh.lduAddr().size(), 0);
        }
        const labelList& fineAgglom =
        (
            nCreatedLevels == startLevel
          ? dummyAgglom
          : agglomeration[nCreatedLevels]
        );
        const label nRegions = max(fineAgglom);

        labelList thisAgglom(mesh.lduAddr().size());

        for (label regioni = 0; regioni < nRegions; regioni++)
        {
            Pout<< "Decomposing region " << regioni << endl;
            CompactListList<label> cellCells;
            const labelList oldToNew
            (
                localCellCells
                (
                    mesh.lduAddr(),
                    fineAgglom,     // marker per cell
                    regioni,        // which marker to keep
                    cellCells
                )
            );

            decomposer.nDomains(4);
            const bool oldParRun = UPstream::parRun(false);
            const labelList cellToProc
            (
                decomposer.decompose
                (
                    cellCells,
                    pointField(cellCells.size(), Zero)
                )
            );
            UPstream::parRun(oldParRun);

            // Combine 
            forAll(oldToNew, celli)
            {
                const label compactCelli = oldToNew[celli];
                if (compactCelli != -1)
                {
                    thisAgglom[celli] = regioni*4 + cellToProc[compactCelli];
                }
            }
        }


        // Store thisAgglom
        agglomeration.append(std::move(thisAgglom));
        nCells_[nCreatedLevels] = max(agglomeration.last())+1;
        restrictAddressing_.set
        (
            nCreatedLevels,
            new labelField(agglomeration.last())
        );

        // Create coarse mesh. Note cannot use agglomerateLduAddressing
        // since going wrong way - it uses fine-mesh

        // Convert current mesh and agglomeration into lduPrimitiveMesh
        {
            const auto& agglom = agglomeration.last();
            const auto& lower = mesh.lduAddr().lowerAddr();
            const auto& upper = mesh.lduAddr().upperAddr();

            labelList faceMap(lower.size(), -1);
            labelList coarseLower(lower.size());
            labelList coarseUpper(lower.size());
            label nCoarseFaces = 0;
            forAll(lower, facei)
            {
                if (faceMap[facei] == -1)
                {
                    const label l = agglom[lower[facei]];
                    const label u = agglom[upper[facei]];

                    if (l < u)
                    {
                        const label coarseFacei = nCoarseFaces++;
                        faceMap[facei] = coarseFacei;
                        coarseLower[coarseFacei] = l;
                        coarseUpper[coarseFacei] = u;
                    }
                    else if (u < l)
                    {
                        const label coarseFacei = nCoarseFaces++;
                        faceMap[facei] = coarseFacei;
                        coarseLower[coarseFacei] = u;
                        coarseUpper[coarseFacei] = l;
                    }
                }
            }
            coarseLower.setSize(nCoarseFaces);
            coarseUpper.setSize(nCoarseFaces);

            Pout<< "Going from" << nl
                << "    nCells:" << mesh.lduAddr().size()
                << " nFaces:" << mesh.lduAddr().lowerAddr().size() << nl
                << "to" << nl
                << "    nCells:" << nCells_[nCreatedLevels]
                << " nFaces:" << coarseLower.size() << endl;

            levelMeshPtr.reset
            (
                new lduPrimitiveMesh
                (
                    nCells_[nCreatedLevels],
                    coarseLower,
                    coarseUpper,
                    mesh.comm(),
                    true
                )
            );
        }
//XXXX
        nCreatedLevels++;
    }


    // Reverse agglomeration table and re-do agglomerateLduAddressing

}
//XXXXXX
// ************************************************************************* //
