/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::masterCoarsestGAMGProcAgglomeration2

Description
    Processor agglomeration of GAMGAgglomerations.

    - by default agglomerates onto the master processor
    - optionally have multiple masters through the nProcessorsPerMaster
      parameter. e.g.

        p
        {
            solver                  GAMG;
            tolerance               1e-06;
            relTol                  0.1;
            smoother                GaussSeidel;

            nCellsInCoarsestLevel   10;

            processorAgglomerator   masterCoarsest2;
            // Groups of 32 cores get combined so each 'master' gets
            // (roughly) 32*nCellsInCoarsestLevel cells.
            nProcessorsPerMaster    32;
        }

SourceFiles
    masterCoarsestGAMGProcAgglomeration2.C

\*---------------------------------------------------------------------------*/

#ifndef masterCoarsestGAMGProcAgglomeration2_H
#define masterCoarsestGAMGProcAgglomeration2_H

#include "GAMGProcAgglomeration.H"
#include "DynamicList.H"
#include "scalarField.H"
#include "labelField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class GAMGAgglomeration;
class lduPrimitiveMesh;
class lduAddressing;

/*---------------------------------------------------------------------------*\
              Class masterCoarsestGAMGProcAgglomeration2 Declaration
\*---------------------------------------------------------------------------*/

class masterCoarsestGAMGProcAgglomeration2
:
    public GAMGProcAgglomeration
{
    // Private data

        const label nProcessorsPerMaster_;

        DynamicList<label> comms_;


    // Private Member Functions

        //- Return (on master) all single-cell meshes collected. single-cell
        //  meshes are just one cell with all proc faces intact.
        autoPtr<lduPrimitiveMesh> singleCellMesh
        (
            const label singleCellMeshComm,
            const lduMesh& mesh,
            scalarField& faceWeights
        ) const;

        tmp<labelField> agglomerate
        (
            label& nCoarseCells,
            const lduAddressing& fineMatrixAddressing,
            const scalarField& faceWeights
        ) const;

        //- Construct processor agglomeration: for every processor the
        //  coarse processor-cluster it agglomerates onto
        tmp<labelField> processorAgglomeration(const lduMesh&) const;

        //- No copy construct
        masterCoarsestGAMGProcAgglomeration2
        (
            const masterCoarsestGAMGProcAgglomeration2&
        ) = delete;

        //- No copy assignment
        void operator=(const masterCoarsestGAMGProcAgglomeration2&) = delete;


public:

    //- Runtime type information
    TypeName("masterCoarsest2");


    // Constructors

        //- Construct given agglomerator and controls
        masterCoarsestGAMGProcAgglomeration2
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    //- Destructor
    virtual ~masterCoarsestGAMGProcAgglomeration2();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
