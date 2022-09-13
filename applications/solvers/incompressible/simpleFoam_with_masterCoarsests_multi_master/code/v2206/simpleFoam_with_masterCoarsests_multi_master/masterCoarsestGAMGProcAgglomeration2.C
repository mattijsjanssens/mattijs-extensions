/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "masterCoarsestGAMGProcAgglomeration2.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"
#include "processorLduInterface.H"
#include "processorGAMGInterface.H"
#include "pairGAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(masterCoarsestGAMGProcAgglomeration2, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        masterCoarsestGAMGProcAgglomeration2,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::lduPrimitiveMesh>
Foam::masterCoarsestGAMGProcAgglomeration2::singleCellMesh
(
    const label singleCellMeshComm,
    const lduMesh& mesh,
    scalarField& faceWeights
) const
{
    // Count number of faces per processor
    List<Map<label>> procFaces(UPstream::nProcs(mesh.comm()));
    Map<label>& myNeighbours = procFaces[UPstream::myProcNo(mesh.comm())];

    {
        const lduInterfacePtrsList interfaces(mesh.interfaces());
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const processorLduInterface& pp =
                    refCast<const processorLduInterface>
                    (
                        interfaces[intI]
                    );
                label size = interfaces[intI].faceCells().size();
                myNeighbours.insert(pp.neighbProcNo(), size);
            }
        }
    }

    Pstream::allGatherList(procFaces, Pstream::msgType(), mesh.comm());

    autoPtr<lduPrimitiveMesh> singleCellMeshPtr;

    if (Pstream::master(mesh.comm()))
    {
        // I am master
        label nCells = Pstream::nProcs(mesh.comm());

        DynamicList<label> l(3*nCells);
        DynamicList<label> u(3*nCells);
        DynamicList<scalar> weight(3*nCells);

        DynamicList<label> nbrs;
        DynamicList<scalar> weights;

        forAll(procFaces, proci)
        {
            const Map<label>& neighbours = procFaces[proci];

            // Add all the higher processors
            nbrs.clear();
            weights.clear();
            forAllConstIters(neighbours, iter)
            {
                if (iter.key() > proci)
                {
                    nbrs.append(iter.key());
                    weights.append(iter());
                }
                sort(nbrs);
                forAll(nbrs, i)
                {
                    l.append(proci);
                    u.append(nbrs[i]);
                    weight.append(weights[i]);
                }
            }
        }

        faceWeights.transfer(weight);

        PtrList<const lduInterface> primitiveInterfaces(0);
        const lduSchedule ps(0);

        singleCellMeshPtr.reset
        (
            new lduPrimitiveMesh
            (
                nCells,
                l,
                u,
                primitiveInterfaces,
                ps,
                singleCellMeshComm
            )
        );
    }
    return singleCellMeshPtr;
}


Foam::tmp<Foam::labelField>
Foam::masterCoarsestGAMGProcAgglomeration2::agglomerate
(
    label& nCoarseCells,
    const lduAddressing& fineMatrixAddressing,
    const scalarField& faceWeights
) const
{
    // Called on master processor only:
    // - agglomerate (multiple steps?)
    // - return number of coarse cells
    // - return map from input to output cell

    return pairGAMGAgglomeration::agglomerate
    (
        nCoarseCells,
        fineMatrixAddressing,
        faceWeights
    );
}


Foam::tmp<Foam::labelField>
Foam::masterCoarsestGAMGProcAgglomeration2::processorAgglomeration
(
    const lduMesh& mesh
) const
{
    label singleCellMeshComm = UPstream::allocateCommunicator
    (
        mesh.comm(),
        labelList(1, Zero)   // only processor 0
    );

    // Return processor-connectivity as a mesh (on master of communicator)
    scalarField faceWeights;
    autoPtr<lduPrimitiveMesh> singleCellMeshPtr
    (
        singleCellMesh
        (
            singleCellMeshComm,
            mesh,
            faceWeights
        )
    );

    tmp<labelField> tfineToCoarse(new labelField(0));
    labelField& fineToCoarse = tfineToCoarse.ref();

    if (singleCellMeshPtr)
    {
        // On master call the agglomerator
        //const lduPrimitiveMesh& singleCellMesh = *singleCellMeshPtr;

        label nCoarseProcs;
        //fineToCoarse = pairGAMGAgglomeration::agglomerate
        fineToCoarse = masterCoarsestGAMGProcAgglomeration2::agglomerate
        (
            nCoarseProcs,
            singleCellMeshPtr(),
            faceWeights
        );

        labelList coarseToMaster(nCoarseProcs, labelMax);
        forAll(fineToCoarse, celli)
        {
            label coarseI = fineToCoarse[celli];
            coarseToMaster[coarseI] = min(coarseToMaster[coarseI], celli);
        }

        // Sort according to master and redo restriction
        labelList newToOld(sortedOrder(coarseToMaster));
        labelList oldToNew(invert(newToOld.size(), newToOld));

        fineToCoarse = labelUIndList(oldToNew, fineToCoarse)();
    }

    Pstream::scatter(fineToCoarse, Pstream::msgType(), mesh.comm());
    UPstream::freeCommunicator(singleCellMeshComm);

    return tfineToCoarse;
}


//XXXXXX
//XXXXXX

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterCoarsestGAMGProcAgglomeration2::masterCoarsestGAMGProcAgglomeration2
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    nProcessorsPerMaster_
    (
        controlDict.getOrDefault<label>("nProcessorsPerMaster", 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterCoarsestGAMGProcAgglomeration2::
~masterCoarsestGAMGProcAgglomeration2()
{
    forAllReverse(comms_, i)
    {
        if (comms_[i] != -1)
        {
            UPstream::freeCommunicator(comms_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::masterCoarsestGAMGProcAgglomeration2::agglomerate()
{
    if (debug)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        // Agglomerate one but last level (since also agglomerating
        // restrictAddressing)
        label fineLevelIndex = agglom_.size()-1;

        if (agglom_.hasMeshLevel(fineLevelIndex))
        {
            // Get the fine mesh
            const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex);
            label levelComm = levelMesh.comm();
            label nProcs = UPstream::nProcs(levelComm);

            if (nProcs > 1)
            {
//                //- Construct given mesh and controls
//                pairGAMGAgglomeration agglomerator
//                (
//                    levelMesh,
//                    dictionary()
//                );


                // Processor restriction map: per processor the coarse processor
                labelList procAgglomMap(nProcs);


DebugVar(nProcessorsPerMaster_);

                if (nProcessorsPerMaster_ > 1)
                {
                    forAll(procAgglomMap, fineProci)
                    {
                        procAgglomMap[fineProci] =
                        (
                            fineProci
                          / nProcessorsPerMaster_
                        );
                    }
                    //procAgglomMap = processorAgglomeration(levelMesh);
                }
                else
                {
                    procAgglomMap = Zero;
                }

DebugVar(procAgglomMap);
Pout<< "procAgglomMap:" << flatOutput(procAgglomMap) << endl;

                // Master processor
                labelList masterProcs;
                // Local processors that agglomerate. agglomProcIDs[0] is in
                // masterProc.
                List<label> agglomProcIDs;
                GAMGAgglomeration::calculateRegionMaster
                (
                    levelComm,
                    procAgglomMap,
                    masterProcs,
                    agglomProcIDs
                );

Pout<< "masterProcs:" << flatOutput(masterProcs) << endl;
Pout<< "agglomProcIDs:" << flatOutput(agglomProcIDs) << endl;

                // Allocate a communicator for the processor-agglomerated matrix
                comms_.append
                (
                    UPstream::allocateCommunicator
                    (
                        levelComm,
                        masterProcs
                    )
                );

                // Use processor agglomeration maps to do the actual collecting.
                if (Pstream::myProcNo(levelComm) != -1)
                {
                    GAMGProcAgglomeration::agglomerate
                    (
                        fineLevelIndex,
                        procAgglomMap,
                        masterProcs,
                        agglomProcIDs,
                        comms_.last()
                    );

                    for
                    (
                        label levelI = 0;
                        levelI <= agglom_.size();
                        levelI++
                    )
                    {
                        if (agglom_.hasMeshLevel(levelI))
                        {
                            const lduMesh& fineMesh = agglom_.meshLevel(levelI);
                            const auto& addr = fineMesh.lduAddr();
                            Pout<< "level:" << levelI
                                << " size:" << addr.size() << endl;


                            const scalarField weights
                            (
                                addr.lowerAddr().size(),
                                1.0
                            );
                            Pout<< "weights:" << weights.size() << endl;
                            
                            dynamic_cast<pairGAMGAgglomeration&>(agglom_)
                            .agglomerate(levelI, weights);
                            
                            Pout<< "**DONE weights:" << weights.size()
                                << endl;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Print a bit
    if (debug)
    {
        Pout<< nl << "Agglomerated mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    return true;
}


// ************************************************************************* //
