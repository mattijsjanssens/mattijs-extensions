/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::kahipGAMGAgglomeration

Description
    Agglomerate using the kahip algorithm.

SourceFiles
    kahipGAMGAgglomeration.C
    kahipGAMGAgglomerate.C

\*---------------------------------------------------------------------------*/

#ifndef kahipGAMGAgglomeration_H
#define kahipGAMGAgglomeration_H

#include "fvMesh.H"
#include "GAMGAgglomeration.H"

// Kahip
#include <stdio.h>

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "timer.h"
#include "coarsening/clustering/size_constraint_label_propagation.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class kahipGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class kahipGAMGAgglomeration
:
    public GAMGAgglomeration
{
    // Private data

        const fvMesh& fvMesh_;


    // Private Member Functions

        void swap
        (
            const lduInterfacePtrsList& interfaces,
            const labelUList& cellValues,
            PtrList<labelList>& nbrValues
        ) const;

        void getNbrAgglom
        (
            const lduAddressing& addr,
            const lduInterfacePtrsList& interfaces,
            const PtrList<labelList>& nbrGlobalAgglom,
            labelList& cellToNbrAgglom
        ) const;

        void detectSharedFaces
        (
            const lduMesh& mesh,
            const labelList& value,
            labelHashSet& sharedFaces
        ) const;


        //- Construct the CSR format addressing
        void makeCompactCellFaceAddressingAndFaceWeights
        (
            const lduAddressing& fineAddressing,
            List<int>& cellCells,
            List<int>& cellCellOffsets,
            const scalarField& magSi,
            List<scalar>& faceWeights
        );

        //- Calculate and return agglomeration
        tmp<labelField> agglomerate
        (
            label& nCoarseCells,
            const label minSize,
            const label maxSize,
            const lduAddressing& fineAddressing,
            const scalarField& V,
            const scalarField& magSf,
            const scalarField& magSb
        );


        //- No copy construct
        kahipGAMGAgglomeration(const kahipGAMGAgglomeration&) = delete;

        //- No copy assignment
        void operator=(const kahipGAMGAgglomeration&) = delete;


public:

    //- Runtime type information
    TypeName("kahip");


    // Constructors

        //- Construct given mesh and controls
        kahipGAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
