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

Class
    Foam::decompositionGAMGAgglomeration

Description
    Agglomerate using the decomposition algorithm

SourceFiles
    decompositionGAMGAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef decompositionGAMGAgglomeration_H
#define decompositionGAMGAgglomeration_H

#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class decompositionMethod;

/*---------------------------------------------------------------------------*\
                 Class decompositionGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class decompositionGAMGAgglomeration
:
    public GAMGAgglomeration
{

    mutable autoPtr<decompositionMethod> decomposerPtr_;


    //- Agglomerate single level
    tmp<labelField> agglomerate
    (
        label& nCoarseCells,
        const lduAddressing& addr,
        const scalarField& faceWeights
    ) const;

    label regionSplit
    (
        const lduAddressing& addr,
        const bitSet& isBlockedFace,
        labelList& cellRegion
    ) const;

    static bitSet blockedFaces
    (
        const lduAddressing& addr,
        const labelUList& region        // region per cell
    );


public:

    //- Runtime type information
    TypeName("decomposition");


    // Constructors

        //- Construct given mesh and controls
        decompositionGAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

        //- Construct given mesh and controls
        decompositionGAMGAgglomeration
        (
            const lduMesh& mesh,
            const scalarField& cellVolumes,
            const vectorField& faceAreas,
            const dictionary& controlDict
        );


    // Member Functions

        //- Agglomerate from a starting level. Starting level is usually 0
        //- (initial mesh) but sometimes >0 (restarting after processor
        //- agglomeration)
        virtual void agglomerate
        (
            const label nCellsInCoarsestLevel,
            const label startLevel,
            const scalarField& startFaceWeights,
            const bool doProcessorAgglomerate = true
        );


    // Helpers

        static void localCellCells
        (
            const lduAddressing& addr,
            const bitSet& isBlockedFace,
            CompactListList<label>& cellCells
        );
        static void localCellCells
        (
            const lduAddressing& addr,
            const labelList& regions,   // marker
            const label regioni,        // which marker to keep
            CompactListList<label>& cellCells
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
