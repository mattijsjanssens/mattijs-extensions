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

Foam::tmp<Foam::labelField> Foam::decompositionGAMGAgglomeration::agglomerate
(
    label& nCoarseCells,
    const lduAddressing& addr,
    const scalarField& faceWeights
)
{
    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

DebugVar(own);
DebugVar(nbr);

    const globalIndex globalCells
    (
        addr.size(),
        UPstream::worldComm,
        false           // local only
    );


    nCoarseCells = 2;

    // Construct globalCellCells for calling decomposition method.
    CompactListList<label> globalCellCells;
    {
        labelList nNbrs(addr.size(), Zero);
        forAll(nbr, facei)
        {
            nNbrs[nbr[facei]]++;
            nNbrs[own[facei]]++;
        }
        globalCellCells.resize_nocopy(nNbrs);

        auto& values = globalCellCells.values();
        auto& offsets = globalCellCells.offsets();

        nNbrs = Zero;
        forAll(nbr, facei)
        {
            const label n = nbr[facei];
            const label o = own[facei];
            values[offsets[n]+(nNbrs[n]++)] = globalCells.toGlobal(o);
            values[offsets[o]+(nNbrs[o]++)] = globalCells.toGlobal(n);
        }
    }


    DebugVar(globalCellCells);


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
DebugVar(nCellsInCoarsestLevel);
DebugVar(startLevel);
DebugVar(startFaceWeights);


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

    //while (nCreatedLevels < maxLevels_ - 1)
    while (nCreatedLevels < startLevel+1)
    {
        DebugVar(nCreatedLevels);

        if (!hasMeshLevel(nCreatedLevels))
        {
            FatalErrorInFunction<< "No mesh at nCreatedLevels:"
                << nCreatedLevels
                << exit(FatalError);
        }

        const auto& fineMesh = meshLevel(nCreatedLevels);


        label nCoarseCells = -1;

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
