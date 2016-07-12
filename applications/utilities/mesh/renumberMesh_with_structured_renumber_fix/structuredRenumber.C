/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "structuredRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "topoDistanceData.H"
#include "fvMeshSubset.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(structuredRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        structuredRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structuredRenumber::structuredRenumber
(
    const dictionary& renumberDict
)
:
    renumberMethod(renumberDict),
    methodDict_(renumberDict.subDict(typeName + "Coeffs")),
    patches_(methodDict_.lookup("patches")),
    //nLayers_(readLabel(methodDict_.lookup("nLayers"))),
    depthFirst_(methodDict_.lookup("depthFirst")),
    method_(renumberMethod::New(methodDict_)),
    reverse_(methodDict_.lookup("reverse"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::structuredRenumber::layerLess::operator()
(
    const label a,
    const label b
)
{
    const topoDistanceData& ta = values_[a];
    const topoDistanceData& tb = values_[b];

    int dummy;

    if (ta.valid(dummy))
    {
        if (tb.valid(dummy))
        {
            if (depthFirst_)
            {
                if (ta.data() < tb.data())
                {
                    // Sort column first
                    return true;
                }
                else if (ta.data() > tb.data())
                {
                    return false;
                }
                else
                {
                    // Same column. Sort according to layer
                    return ta.distance() < tb.distance();
                }
            }
            else
            {
                if (ta.distance() < tb.distance())
                {
                    return true;
                }
                else if (ta.distance() > tb.distance())
                {
                    return false;
                }
                else
                {
                    // Same layer; sort according to column
                    return ta.data() < tb.data();
                }
            }
        }
        else
        {
            return true;
        }
    }
    else
    {
        if (tb.valid(dummy))
        {
            return false;
        }
        else
        {
            // Both not valid; fall back to cell indices for sorting
            return a < b;
        }
    }
}


Foam::labelList Foam::structuredRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Number of points " << points.size()
            << " should equal the number of cells " << mesh.nCells()
            << exit(FatalError);
    }

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const labelHashSet patchIDs(pbm.patchSet(patches_));

    label nFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nFaces += pbm[iter.key()].size();
    }


    // Extract a submesh.
    labelHashSet patchCells(2*nFaces);
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const labelUList& fc = pbm[iter.key()].faceCells();
        forAll(fc, i)
        {
            patchCells.insert(fc[i]);
        }
    }

    label nTotalSeeds = returnReduce(patchCells.size(), sumOp<label>());

    label nTotalCells = mesh.globalData().nTotalCells();
    const label nLayers = nTotalCells/nTotalSeeds;

    Info<< type() << " : seeding " << nTotalSeeds
        << " cells on (estimated) " << nLayers << " layers" << nl
        << endl;


    // Avoid subsetMesh, FaceCellWave going through proc boundaries
    // Note: not too sure that FaceCellWave shouldn't since e.g. a processor
    // might have zero patch faces on it but be quite close to the patch.
    // However effect of renumbering is always local since processor boundaries
    // are explicit.

//    bool oldParRun = Pstream::parRun();
//    Pstream::parRun() = false;


    // Work array. Used here to temporarily store the original-to-ordered
    // index. Later on used to store the ordered-to-original.
    labelList orderedToOld(mesh.nCells(), -1);

    // Subset the layer of cells next to the patch
    {
        fvMeshSubset subsetter(dynamic_cast<const fvMesh&>(mesh));
        subsetter.setLargeCellSubset(patchCells);
        const fvMesh& subMesh = subsetter.subMesh();

        pointField subPoints(points, subsetter.cellMap());

        // Locally renumber the layer of cells
        labelList subOrder(method_().renumber(subMesh, subPoints));

        labelList subOrigToOrdered(invert(subOrder.size(), subOrder));

        globalIndex globalSubCells(subOrder.size());

        // Transfer to final decomposition and convert into global numbering
        forAll(subOrder, i)
        {
            orderedToOld[subsetter.cellMap()[i]] =
                globalSubCells.toGlobal(subOrigToOrdered[i]);
        }
    }


    // Walk sub-ordering (=column index) out.
    labelList patchFaces(nFaces);
    List<topoDistanceData> patchData(nFaces);
    nFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const polyPatch& pp = pbm[iter.key()];
        const labelUList& fc = pp.faceCells();
        forAll(fc, i)
        {
            patchFaces[nFaces] = pp.start()+i;
            patchData[nFaces] = topoDistanceData
            (
                orderedToOld[fc[i]],// passive data: global column
                0                   // distance: layer
            );
            nFaces++;
        }
    }

    // Field on cells and faces.
    List<topoDistanceData> cellData(mesh.nCells());
    List<topoDistanceData> faceData(mesh.nFaces());

    // Propagate information inwards
    FaceCellWave<topoDistanceData> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        faceData,
        cellData,
        nTotalCells+1
    );


//     Pstream::parRun() = oldParRun;

    // We now use the global column and layer in column to determine a sorting
    // order. Combine column and layer into single index
    // (alternatively we could use a specialised comparison operator)

//     // Determine max layer (max column is already known = nTotalSeeds)
//     label maxLayer = 0;
//     forAll(cellData, celli)
//     {
//         if (cellData[celli].valid(deltaCalc.data()))
//         {
//             maxLayer = max(maxLayer, cellData[celli].distance());
//         }
//     }
// 
//     // Determine max layers
// 
// 
//     // And extract
//     label nUnvisited = 0;
//     labelList cellDistance(mesh.nCells(), 0);
//     forAll(cellData, celli)
//     {
//         if (!cellData[celli].valid(deltaCalc.data()))
//         {
//             nUnvisited++;
//         }
//         else
//         {
//             label columnI = cellData[celli].data();
//             label layerI = cellData[celli].distance();
// 
//             if (depthFirst_)
//             {
//                 cellDistance[celli] = nTotalSeeds*columnI+layerI;
//             }
//             else
//             {
//                 cellDistance[celli] = columnI+nTotalSeeds*layerI;
//             }
//         }
//     }

    Info<< type() << " : did not visit "
        << deltaCalc.getUnsetCells()
        << " cells out of " << nTotalCells
        << "; keeping these in original order" << endl;


    //labelList orderedToOld;
    sortedOrder(cellData, orderedToOld, layerLess(depthFirst_, cellData));

// 
//     // Note that distance is distance from face so starts at 1.
//     bool haveWarned = false;
//     forAll(orderedToOld, celli)
//     {
//         if (!cellData[celli].valid(deltaCalc.data()))
//         {
//             if (!haveWarned)
//             {
//                 WarningInFunction
//                     << "Did not visit some cells, e.g. cell " << celli
//                     << " at " << mesh.cellCentres()[celli] << endl
//                     << "Assigning these cells to domain 0." << endl;
//                 haveWarned = true;
//             }
//             orderedToOld[celli] = 0;
//         }
//         else
//         {
//             label layerI = cellData[celli].distance();
//             if (depthFirst_)
//             {
//                 orderedToOld[nLayers*cellData[celli].data()+layerI] = celli;
//             }
//             else
//             {
//                 orderedToOld[cellData[celli].data()+nLayers*layerI] = celli;
//             }
//         }
//     }

    // Return furthest away cell first
    if (reverse_)
    {
        reverse(orderedToOld);
    }

    return orderedToOld;
}


// ************************************************************************* //
