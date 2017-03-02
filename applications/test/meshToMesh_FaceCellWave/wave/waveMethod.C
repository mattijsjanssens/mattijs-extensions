/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "waveMethod.H"
#include "meshToMeshData.H"
#include "FaceCellWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(waveMethod, 0);
    addToRunTimeSelectionTable(meshToMeshMethod, waveMethod, components);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::waveMethod::calculate
(
    const polyMesh& src,
    const polyMesh& tgt,
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght
) const
{
    DynamicList<label> changedFaces(src.nFaces()/100 + 100);
    DynamicList<meshToMeshData> changedFacesInfo(changedFaces.size());

    List<meshToMeshData> cellData(src.nCells());
    List<meshToMeshData> faceData(src.nFaces());

    meshToMeshData::trackData td(tgt);

    while (true)
    {
        changedFaces.clear();
        changedFacesInfo.clear();

        // Search for starting seed
        forAll(cellData, celli)
        {
            if (!cellData[celli].valid(td))
            {
                const point& cc = src.cellCentres()[celli];
                label tgtCelli = td.tgtMesh_.findCell(cc, polyMesh::CELL_TETS);
                if (tgtCelli != -1)
                {
                    // Insert any face of cell
                    label facei = src.cells()[celli][0];
                    changedFaces.append(facei);
                    changedFacesInfo.append(meshToMeshData(tgtCelli));
                    break;
                }
                else
                {
                    // Register as no correspondence
                    cellData[celli] = meshToMeshData(-1);
                }
            }
        }

        if (changedFaces.empty())
        {
            break;
        }

        FaceCellWave<meshToMeshData, meshToMeshData::trackData> calc
        (
            src,
            changedFaces,
            changedFacesInfo,
            faceData,
            cellData,
            src.globalData().nTotalCells(),   // max iterations
            td
        );
    }


    // Copy into srcToTgt
    srcToTgtAddr.setSize(src.nCells());
    srcToTgtWght.setSize(src.nCells());

    forAll(cellData, celli)
    {
        if (cellData[celli].tgtCelli_ != -1)
        {
            labelList& l = srcToTgtAddr[celli];
            scalarList& w = srcToTgtWght[celli];
            l.setSize(1);
            w.setSize(1);
            l[0] = cellData[celli].tgtCelli_;
            w[0] = src.cellVolumes()[celli];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveMethod::waveMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    meshToMeshMethod(src, tgt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveMethod::~waveMethod()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveMethod::calculate
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght
)
{
    calculate(src_, tgt_, srcToTgtAddr, srcToTgtWght);
    calculate(tgt_, src_, tgtToSrcAddr, tgtToSrcWght);
}


// ************************************************************************* //
