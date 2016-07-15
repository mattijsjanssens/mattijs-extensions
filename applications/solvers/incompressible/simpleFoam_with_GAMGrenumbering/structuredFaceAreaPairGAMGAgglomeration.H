/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

Class
    Foam::structuredFaceAreaPairGAMGAgglomeration

Description
    Agglomerate using the pair algorithm.

SourceFiles
    structuredFaceAreaPairGAMGAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef structuredFaceAreaPairGAMGAgglomeration_H
#define structuredFaceAreaPairGAMGAgglomeration_H

#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class structuredFaceAreaPairGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class structuredFaceAreaPairGAMGAgglomeration
:
    public GAMGAgglomeration
{
    // Private data

        //- Number of levels to merge, 1 = don't merge, 2 = merge pairs etc.
        label mergeLevels_;

        //- Direction of cell loop for the current level
        static bool forward_;

protected:

    // Protected Member Functions

        //- Agglomerate all levels starting from the given face weights
        void agglomerate
        (
            const lduMesh& mesh,
            const scalarField& faceWeights
        );

        //- Calculate and return agglomeration
        static tmp<labelField> structuredAgglomerate
        (
            label& nCoarseCells,
            const lduAddressing& fineMatrixAddressing,
            const scalarField& faceWeights
        );

        //- Calculate and return agglomeration
        static tmp<labelField> agglomerate
        (
            label& nCoarseCells,
            const lduAddressing& fineMatrixAddressing,
            const scalarField& faceWeights
        );

        //- Disallow default bitwise copy construct
        structuredFaceAreaPairGAMGAgglomeration(const structuredFaceAreaPairGAMGAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const structuredFaceAreaPairGAMGAgglomeration&);


public:

    //- Runtime type information
    TypeName("structuredFaceAreaPair");


    // Constructors

        //- Construct given mesh and controls
        structuredFaceAreaPairGAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

//         //- Construct given mesh and controls
//         structuredFaceAreaPairGAMGAgglomeration
//         (
//             const lduMesh& mesh,
//             const scalarField& cellVolumes,
//             const vectorField& faceAreas,
//             const dictionary& controlDict
//         );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
