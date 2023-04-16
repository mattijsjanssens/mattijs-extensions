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

Description
    Agglomerate one level using the kahip algorithm.

\*---------------------------------------------------------------------------*/

#include "kahipGAMGAgglomeration.H"
//#include "fvMesh.H"
#include "OFstream.H"

#include "random_functions.h"
#include "graph_io.h"
#include "configuration.h"
//#include "parse_parameters.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::kahipGAMGAgglomeration::
makeCompactCellFaceAddressingAndFaceWeights
(
    const lduAddressing& fineAddressing,
    List<int>& cellCells,
    List<int>& cellCellOffsets,
    const scalarField& magSi,
    List<scalar>& faceWeights
)
{
    const label nFineCells = fineAddressing.size();
    const label nFineFaces = fineAddressing.upperAddr().size();

    const labelUList& upperAddr = fineAddressing.upperAddr();
    const labelUList& lowerAddr = fineAddressing.lowerAddr();

    // Number of neighbours for each cell
    labelList nNbrs(nFineCells, Zero);

    forAll(upperAddr, facei)
    {
        nNbrs[upperAddr[facei]]++;
    }

    forAll(lowerAddr, facei)
    {
        nNbrs[lowerAddr[facei]]++;
    }

    // Set the sizes of the addressing and faceWeights arrays
    cellCellOffsets.setSize(nFineCells + 1);
    cellCells.setSize(2*nFineFaces);
    faceWeights.setSize(2*nFineFaces);


    cellCellOffsets[0] = 0;
    forAll(nNbrs, celli)
    {
        cellCellOffsets[celli+1] = cellCellOffsets[celli] + nNbrs[celli];
    }

    // reset the whole list to use as counter
    nNbrs = 0;

    forAll(upperAddr, facei)
    {
        label own = upperAddr[facei];
        label nei = lowerAddr[facei];

        label l1 = cellCellOffsets[own] + nNbrs[own]++;
        label l2 = cellCellOffsets[nei] + nNbrs[nei]++;

        cellCells[l1] = nei;
        cellCells[l2] = own;

        faceWeights[l1] = magSi[facei];
        faceWeights[l2] = magSi[facei];
    }
}


Foam::tmp<Foam::labelField> Foam::kahipGAMGAgglomeration::agglomerate
(
    label& nCoarseCells,
    const label minSize,
    const label maxSize,
    const lduAddressing& fineAddressing,
    const scalarField& V,
    const scalarField& magSf,
    const scalarField& magSb
)
{
    const label nFineCells = fineAddressing.size();

    // Compact addressing for cellCells
    List<int> cellCells;
    List<int> cellCellOffsets;

    // Face weights = face areas of the internal faces
    List<scalar> faceWeights;

    // Create the compact addressing for cellCells and faceWeights
    makeCompactCellFaceAddressingAndFaceWeights
    (
        fineAddressing,
        cellCells,
        cellCellOffsets,
        magSf,
        faceWeights
    );

    if (debug)
    {
        const label nFineFaces = fineAddressing.upperAddr().size();
        const fileName graph_filename("graph_" + Foam::name(nFineCells) + ".grf");

        {
            OFstream os(graph_filename);
            Pout<< "Dumping Metis graph file for nCells "
                << nFineCells << " to " << os.name() << endl;


            os << "% Metis graph format" << nl
                << nFineCells << ' ' << nFineFaces << nl;

            for (label celli = 0; celli < nFineCells; celli++)
            {
                const label start = cellCellOffsets[celli];
                const label end = cellCellOffsets[celli+1];

                // First nbr
                os << cellCells[start]+1;

                // Rest of nbrs
                const SubList<label> nbrs(cellCells, end-start-1, start+1);

                for (const label nbr : nbrs)
                {
                    os  << ' ' << nbr+1;
                }
                os << nl;
            }
        }

if (false)
        {
DebugVar("**READ**");
            graph_access G;
            graph_io::readGraphWeighted(G, graph_filename);

            PartitionConfig partition_config;

            configuration cfg;
            cfg.standard(partition_config);
            cfg.fast_separator(partition_config);
            partition_config.cluster_upperbound =
                std::numeric_limits< NodeWeight >::max()/2;
//
//            bool is_graph_weighted = false;
//            bool suppress_output   = false;
//            bool recursive         = false;
//
//            int argn = 0;
//            char **argv = nullptr;
//            int ret_code = parse_parameters(argn, argv, 
//                                            partition_config, 
//                                            graph_filename, 
//                                            is_graph_weighted, 
//                                            suppress_output, recursive); 



            std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
            if( partition_config.cluster_upperbound == std::numeric_limits< NodeWeight >::max()/2 ) {
                    std::cout <<  "no size-constrained specified" << std::endl;
            } else {
                    std::cout <<  "size-constrained set to " <<  partition_config.cluster_upperbound << std::endl;
            }

            partition_config.upper_bound_partition = partition_config.cluster_upperbound+1;
            partition_config.cluster_coarsening_factor = 1;
            partition_config.k = 1;
            srand(partition_config.seed);
            random_functions::setSeed(partition_config.seed);

            // ***************************** perform clustering ***************************************       
            NodeID no_blocks = 0;
            std::vector< NodeID > cluster_id(G.number_of_nodes());
            size_constraint_label_propagation sclp;
            sclp.label_propagation( partition_config, G, cluster_id, no_blocks);
            Pout<< "no_blocks:" << no_blocks << endl;

            // output some information about the partition that we have computed 
            quality_metrics qm;
            forall_nodes(G, node) {
                    G.setPartitionIndex(node, cluster_id[node]);
            } endfor

            G.set_partition_count(no_blocks);
            std::cout << "number of clusters/blocks  " << no_blocks << std::endl;
            std::cout << "number of edges between clusters " << qm.edge_cut(G)                 << std::endl;

DebugVar("**END OF READ**");
        }
    }



DebugVar("**PROGRAM**");

    PartitionConfig partition_config;
    partition_config.k = 1;        // wanted cluster size
    configuration cfg;
    cfg.standard(partition_config);
    cfg.fast(partition_config);

    partition_config.cluster_upperbound =
        std::numeric_limits< NodeWeight >::max()/2;

    graph_access G;
    //internal_build_graph
    //(
    //    partition_config,
    //    nFineCells,
    //    nullptr,    // vertex weights
    //    cellCellOffsets.begin(),    // offsets
    //    nullptr,    // adjcwgt,
    //    cellCells.begin(),  //adjncy,
    //    G
    //);
    G.build_from_metis(nFineCells, cellCellOffsets.begin(), cellCells.begin()); 
    //G.set_partition_count(partition_config.k); 
    //if(vwgt != NULL) {
    //        forall_nodes(G, node) {
    //                G.setNodeWeight(node, vwgt[node]);
    //        } endfor
    //}
    //
    //if(adjcwgt != NULL) {
    //        forall_edges(G, e) {
    //                G.setEdgeWeight(e, adjcwgt[e]);
    //        } endfor 
    //}

    std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and "
        <<  G.number_of_edges() <<  " edges"  << std::endl;


//    partition_config.cluster_upperbound =
//        std::numeric_limits< NodeWeight >::max()/2;

    if( partition_config.cluster_upperbound == std::numeric_limits< NodeWeight >::max()/2 ) {
            std::cout <<  "no size-constrained specified" << std::endl;
    } else {
            std::cout <<  "size-constrained set to " <<  partition_config.cluster_upperbound << std::endl;
    }

    partition_config.upper_bound_partition =
        partition_config.cluster_upperbound+1;
    partition_config.cluster_coarsening_factor = 1;
    partition_config.k = 1;        // wanted cluster size
    //partition_config.seed = 0;
    srand(partition_config.seed);
    random_functions::setSeed(partition_config.seed);

    //configuration cfg;
    //cfg.fast(partition_config);     // use 'fast' algo

    NodeID no_blocks = 0;
    std::vector< NodeID > cluster_id(G.number_of_nodes());
    size_constraint_label_propagation sclp;
    sclp.label_propagation( partition_config, G, cluster_id, no_blocks);

    DebugVar(no_blocks);


    forall_nodes(G, node)
    {
        Pout<< "    node:" << node << " cluster:" << cluster_id[node] << endl;
        G.setPartitionIndex(node, cluster_id[node]);
    }
    endfor
    G.set_partition_count(no_blocks);

    nCoarseCells = no_blocks;
    List<int> finalAgglom(nFineCells);
    forall_nodes(G, node)
    {
        finalAgglom[node] = cluster_id[node];
    }
    endfor

//    {
//        label nNewCoarseCells = 0;
//        labelList newRestrictAddr;
//        bool ok = checkRestriction
//        (
//            newRestrictAddr,
//            nNewCoarseCells,
//            fineAddressing,
//            finalAgglom,
//            nCoarseCells
//        );
//
//        if (!ok)
//        {
//            nCoarseCells = nNewCoarseCells;
//            finalAgglom.transfer(newRestrictAddr);
//        }
//    }

    return tmp<labelField>::New(finalAgglom);
}


// ************************************************************************* //
