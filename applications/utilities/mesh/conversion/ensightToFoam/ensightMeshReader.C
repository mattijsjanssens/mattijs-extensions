/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "ensightMeshReader.H"
//#include "oldCyclicPolyPatch.H"
//#include "emptyPolyPatch.H"
//#include "wallPolyPatch.H"
//#include "symmetryPolyPatch.H"
#include "cellModel.H"
//#include "ListOps.H"
//#include "stringOps.H"
//#include "IFstream.H"
//#include "IOMap.H"
#include "ensightReadFile.H"
#include "OBJstream.H"
#include "matchPoints.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//    // Read and discard to newline
//    static inline void readToNewline(ISstream& is)
//    {
//        char ch = '\n';
//        do
//        {
//            is.get(ch);
//        }
//        while ((is) && ch != '\n');
//    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//void Foam::fileFormats::ensightMeshReader::readAux
//(
//    const objectRegistry& registry
//)
//{
//    boundaryRegion_.readDict(registry);
//    cellTable_.readDict(registry);
//}


// Read points from <.vrt> file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_VERTEX [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <vertexId>  <x>  <y>  <z> [newline]

\*---------------------------------------------------------------------------*/
//Foam::label Foam::fileFormats::ensightMeshReader::readPoints
//(
//    const fileName& inputName,
//    const scalar scaleFactor
//)
//{
//    label nPoints = 0, maxId = 0;
//    token tok;
//
//    // Pass 1:
//    // get # points and maximum vertex label
//    {
//        IFstream is(inputName);
//        readHeader(is, ensightCore::HEADER_VRT);
//
//        scalar x, y, z;
//
//        while (is.read(tok).good() && tok.isLabel())
//        {
//            const label starVertexId = tok.labelToken();
//
//            is >> x >> y >> z;
//
//            maxId = max(maxId, starVertexId);
//            ++nPoints;
//        }
//    }
//
//    if (!nPoints)
//    {
//        FatalErrorInFunction
//            << "No points in file " << inputName << nl
//            << abort(FatalError);
//    }
//
//    Info<< "Number of points  = " << nPoints << endl;
//
//    // Set sizes and reset to invalid values
//
//    points_.setSize(nPoints);
//    mapToFoamPointId_.setSize(maxId+1);
//
//    //- Original Point number for a given vertex
//    // might need again in the future
//    ////     labelList origPointId(nPoints);
//    ////     origPointId = -1;
//
//    mapToFoamPointId_ = -1;
//
//    // Pass 2:
//    // construct pointList and conversion table
//    // from Star vertex numbers to Foam point labels
//    {
//        IFstream is(inputName);
//        readHeader(is, ensightCore::HEADER_VRT);
//
//        label pointi = 0;
//        while (is.read(tok).good() && tok.isLabel())
//        {
//            const label starVertexId = tok.labelToken();
//
//            is  >> points_[pointi].x()
//                >> points_[pointi].y()
//                >> points_[pointi].z();
//
//            // might need again in the future
//            ////  origPointId[pointi] = starVertexId;
//            mapToFoamPointId_[starVertexId] = pointi;
//            ++pointi;
//        }
//
//        if (nPoints > pointi)
//        {
//            nPoints = pointi;
//            points_.setSize(nPoints);
//            // might need again in the future
//            //// origPointId.setSize(nPoints);
//        }
//
//        if
//        (
//            scaleFactor > 0
//         && (scaleFactor > 1.0 + SMALL || scaleFactor < 1.0 - SMALL)
//        )
//        {
//            points_ *= scaleFactor;
//        }
//    }
//
//    return maxId;
//}


// Read cells from <.cel> file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_CELL [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <cellId>  <shapeId>  <nLabels>  <cellTableId>  <typeId> [newline]
  <cellId>  <int1> .. <int8>
  <cellId>  <int9> .. <int16>

 with shapeId:
 *   1 = point
 *   2 = line
 *   3 = shell
 *  11 = hexa
 *  12 = prism
 *  13 = tetra
 *  14 = pyramid
 * 255 = polyhedron

 with typeId
 *   1 = fluid
 *   2 = solid
 *   3 = baffle
 *   4 = shell
 *   5 = line
 *   6 = point

For primitive cell shapes, the number of vertices will never exceed 8 (hexa)
and corresponds to <nLabels>.
For polyhedral, <nLabels> includes an index table comprising beg/end pairs
for each cell face.

Strictly speaking, we only need the cellModeller for adding boundaries.
\*---------------------------------------------------------------------------*/

//void Foam::fileFormats::ensightMeshReader::readCells(const fileName& inputName)
//{
//    label nFluids = 0, nSolids = 0, nBaffles = 0, nShells = 0;
//    label maxId = 0;
//    token tok;
//
//    bool unknownVertices = false;
//
//
//    // Pass 1:
//    // count nFluids, nSolids, nBaffle, nShell and maxId
//    // also see if polyhedral cells were used
//    {
//        IFstream is(inputName);
//        readHeader(is, ensightCore::HEADER_CEL);
//
//        label shapeId, nLabels, cellTableId, typeId;
//
//        while (is.read(tok).good() && tok.isLabel())
//        {
//            const label starCellId = tok.labelToken();
//
//            is  >> shapeId
//                >> nLabels
//                >> cellTableId
//                >> typeId;
//
//            // Skip the rest of the line
//            readToNewline(is);
//
//            // Max 8 indices per line
//            while (nLabels > 0)
//            {
//                readToNewline(is);
//                nLabels -= 8;
//            }
//
//            if (typeId == ensightCore::starcdFluidType)
//            {
//                ++nFluids;
//                maxId = max(maxId, starCellId);
//
//                if (!cellTable_.found(cellTableId))
//                {
//                    cellTable_.setName(cellTableId);
//                    cellTable_.setMaterial(cellTableId, "fluid");
//                }
//            }
//            else if (typeId == ensightCore::starcdSolidType)
//            {
//                ++nSolids;
//                if (keepSolids_)
//                {
//                    maxId = max(maxId, starCellId);
//                }
//
//                if (!cellTable_.found(cellTableId))
//                {
//                    cellTable_.setName(cellTableId);
//                    cellTable_.setMaterial(cellTableId, "solid");
//                }
//            }
//            else if (typeId == ensightCore::starcdBaffleType)
//            {
//                // baffles have no cellTable entry
//                ++nBaffles;
//                maxId = max(maxId, starCellId);
//            }
//            else if (typeId == ensightCore::starcdShellType)
//            {
//                ++nShells;
//                if (!cellTable_.found(cellTableId))
//                {
//                    cellTable_.setName(cellTableId);
//                    cellTable_.setMaterial(cellTableId, "shell");
//                }
//            }
//        }
//    }
//
//    const label nCells = nFluids + (keepSolids_ ? nSolids : 0);
//
//    Info<< "Number of fluids  = " << nFluids << nl
//        << "Number of baffles = " << nBaffles << nl
//        << "Number of solids  = " << nSolids
//        << (keepSolids_ ? " (treat as fluid)" : " (ignored)") << nl
//        << "Number of shells  = " << nShells << " (ignored)" << nl;
//
//    if (!nCells)
//    {
//        OSstream& err = FatalErrorInFunction;
//
//        err << "No cells in file " << inputName << nl;
//
//        if (nShells)
//        {
//            err << "Consists of shells only (typeId=4)." << nl;
//        }
//
//        err << nl
//            << abort(FatalError);
//    }
//
//    cellFaces_.setSize(nCells);
//    cellShapes_.setSize(nCells);
//    cellTableId_.setSize(nCells);
//
//    // information for the interfaces
//    baffleFaces_.setSize(nBaffles);
//
//    // extra space for baffles
//    origCellId_.setSize(nCells + nBaffles);
//    mapToFoamCellId_.setSize(maxId+1);
//    mapToFoamCellId_ = -1;
//
//
//    // avoid undefined shapes for polyhedra
//    cellShape genericShape
//    (
//        cellModel::ref(cellModel::UNKNOWN), labelList()
//    );
//
//    // Pass 2:
//    // construct cellFaces_ and possibly cellShapes_
//    {
//        IFstream is(inputName);
//        readHeader(is, ensightCore::HEADER_CEL);
//
//        labelList starLabels(64);
//        label ignoredLabel, shapeId, nLabels, cellTableId, typeId;
//
//        label celli = 0, bafflei = 0;
//
//        while (is.read(tok).good() && tok.isLabel())
//        {
//            const label starCellId = tok.labelToken();
//
//            is  >> shapeId
//                >> nLabels
//                >> cellTableId
//                >> typeId;
//
//            if (nLabels > starLabels.size())
//            {
//                starLabels.setSize(nLabels);
//            }
//            starLabels = -1;
//
//            // Read indices - max 8 per line
//            for (label i = 0; i < nLabels; ++i)
//            {
//                if ((i % 8) == 0)
//                {
//                    is >> ignoredLabel; // Skip cellId for continuation lines
//                }
//                is >> starLabels[i];
//            }
//
//            // Skip solid cells
//            if
//            (
//                typeId == ensightCore::starcdSolidType
//             && !keepSolids_
//            )
//            {
//                continue;
//            }
//
//            // Determine the OpenFOAM cell shape
//            const cellModel* curModelPtr = nullptr;
//
//            // fluid/solid cells
//            switch (shapeId)
//            {
//                case ensightCore::starcdHex:
//                    curModelPtr = cellModel::ptr(cellModel::HEX);
//                    break;
//                case ensightCore::starcdPrism:
//                    curModelPtr = cellModel::ptr(cellModel::PRISM);
//                    break;
//                case ensightCore::starcdTet:
//                    curModelPtr = cellModel::ptr(cellModel::TET);
//                    break;
//                case ensightCore::starcdPyr:
//                    curModelPtr = cellModel::ptr(cellModel::PYR);
//                    break;
//            }
//
//            if (curModelPtr)
//            {
//                // primitive cell - use shapes
//
//                // convert orig vertex Id to point label
//                bool isBad = false;
//                for (label i=0; i < nLabels; ++i)
//                {
//                    label pointId = mapToFoamPointId_[starLabels[i]];
//                    if (pointId < 0)
//                    {
//                        Info<< "Cells inconsistent with vertex file. "
//                            << "Star vertex " << starLabels[i]
//                            << " does not exist" << endl;
//                        isBad = true;
//                        unknownVertices = true;
//                    }
//                    starLabels[i] = pointId;
//                }
//
//                if (isBad)
//                {
//                    continue;
//                }
//
//                // record original cell number and lookup
//                origCellId_[celli] = starCellId;
//                mapToFoamCellId_[starCellId] = celli;
//
//                cellTableId_[celli] = cellTableId;
//                cellShapes_[celli] = cellShape
//                (
//                    *curModelPtr,
//                    SubList<label>(starLabels, nLabels)
//                );
//
//                cellFaces_[celli] = cellShapes_[celli].faces();
//                ++celli;
//            }
//            else if (shapeId == ensightCore::starcdPoly)
//            {
//                // polyhedral cell
//                label nFaces = starLabels[0] - 1;
//
//                // convert orig vertex id to point label
//                // start with offset (skip the index table)
//                bool isBad = false;
//                for (label i=starLabels[0]; i < nLabels; ++i)
//                {
//                    label pointId = mapToFoamPointId_[starLabels[i]];
//                    if (pointId < 0)
//                    {
//                        Info<< "Cells inconsistent with vertex file. "
//                            << "Star vertex " << starLabels[i]
//                            << " does not exist" << endl;
//                        isBad = true;
//                        unknownVertices = true;
//                    }
//                    starLabels[i] = pointId;
//                }
//
//                if (isBad)
//                {
//                    continue;
//                }
//
//                // traverse beg/end indices
//                faceList faces(nFaces);
//                label facei = 0;
//                for (label i=0; i < nFaces; ++i)
//                {
//                    label beg = starLabels[i];
//                    label n   = starLabels[i+1] - beg;
//
//                    face f
//                    (
//                        SubList<label>(starLabels, n, beg)
//                    );
//
//                    f.collapse();
//
//                    // valid faces only
//                    if (f.size() >= 3)
//                    {
//                        faces[facei++] = f;
//                    }
//                }
//
//                if (nFaces > facei)
//                {
//                    Info<< "star cell " << starCellId << " has "
//                        << (nFaces - facei)
//                        << " empty faces - could cause boundary "
//                        << "addressing problems"
//                        << endl;
//
//                    nFaces = facei;
//                    faces.setSize(nFaces);
//                }
//
//                if (nFaces < 4)
//                {
//                    FatalErrorInFunction
//                        << "star cell " << starCellId << " has " << nFaces
//                        << abort(FatalError);
//                }
//
//                // record original cell number and lookup
//                origCellId_[celli] = starCellId;
//                mapToFoamCellId_[starCellId] = celli;
//
//                cellTableId_[celli] = cellTableId;
//                cellShapes_[celli]  = genericShape;
//                cellFaces_[celli]   = faces;
//                ++celli;
//            }
//            else if (typeId == ensightCore::starcdBaffleType)
//            {
//                // baffles
//
//                // convert orig vertex id to point label
//                bool isBad = false;
//                for (label i=0; i < nLabels; ++i)
//                {
//                    label pointId = mapToFoamPointId_[starLabels[i]];
//                    if (pointId < 0)
//                    {
//                        Info<< "Baffles inconsistent with vertex file. "
//                            << "Star vertex " << starLabels[i]
//                            << " does not exist" << endl;
//                        isBad = true;
//                        unknownVertices = true;
//                    }
//                    starLabels[i] = pointId;
//                }
//
//                if (isBad)
//                {
//                    continue;
//                }
//
//
//                face f
//                (
//                    SubList<label>(starLabels, nLabels)
//                );
//
//                f.collapse();
//
//                // valid faces only
//                if (f.size() >= 3)
//                {
//                    baffleFaces_[bafflei] = f;
//                    // insert lookup addressing in normal list
//                    mapToFoamCellId_[starCellId]  = nCells + bafflei;
//                    origCellId_[nCells + bafflei] = starCellId;
//                    ++bafflei;
//                }
//            }
//        }
//
//        baffleFaces_.setSize(bafflei);
//    }
//
//    if (unknownVertices)
//    {
//        FatalErrorInFunction
//            << "cells with unknown vertices"
//            << abort(FatalError);
//    }
//
//    // truncate lists
//
//#ifdef DEBUG_READING
//    Info<< "CELLS READ" << endl;
//#endif
//
//    // cleanup
//    mapToFoamPointId_.clear();
//}


// Read boundaries from <.bnd> file
//
/*---------------------------------------------------------------------------*\
Line 1:
  PROSTAR_BOUNDARY [newline]

Line 2:
  <version> 0 0 0 0 0 0 0 [newline]

Body:
  <boundId>  <cellId>  <cellFace>  <regionId>  0  <boundaryType> [newline]

where boundaryType is truncated to 4 characters from one of the following:
INLET
PRESSSURE
OUTLET
BAFFLE
etc,
\*---------------------------------------------------------------------------*/

//void Foam::fileFormats::ensightMeshReader::readBoundary
//(
//    const fileName& inputName
//)
//{
//    label nPatches = 0, nFaces = 0, nBafflePatches = 0, maxId = 0;
//    label starCellId, cellFaceId, starRegion, configNumber;
//    token tok;
//    word patchType;
//
//    labelList mapToFoamPatchId(1000, label(-1));
//    labelList nPatchFaces(1000, Zero);
//    labelList origRegion(1000, Zero);
//    patchTypes_.setSize(1000);
//
//    //
//    // Mapping between OpenFOAM and PROSTAR primitives
//    // - needed for face mapping
//    //
//    const Map<label> shapeLookup =
//    {
//        { cellModel::ref(cellModel::HEX).index(), ensightCore::starcdHex },
//        { cellModel::ref(cellModel::PRISM).index(), ensightCore::starcdPrism },
//        { cellModel::ref(cellModel::TET).index(), ensightCore::starcdTet },
//        { cellModel::ref(cellModel::PYR).index(), ensightCore::starcdPyr },
//    };
//
//    // Pass 1:
//    // collect
//    // no. of faces (nFaces), no. of patches (nPatches)
//    // and for each of these patches the number of faces
//    // (nPatchFaces[patchLabel])
//    //
//    // and a conversion table from Star regions to (Foam) patchLabels
//    //
//    // additionally note the no. of baffle patches (nBafflePatches)
//    // so that we sort these to the end of the patch list
//    // - this makes it easier to transfer them to an adjacent patch if reqd
//    {
//        IFstream is(inputName);
//
//        if (is.good())
//        {
//            readHeader(is, ensightCore::HEADER_BND);
//
//            while (is.read(tok).good() && tok.isLabel())
//            {
//                // Ignore boundary id (not needed)
//
//                ++nFaces;
//
//                is  >> starCellId
//                    >> cellFaceId
//                    >> starRegion
//                    >> configNumber
//                    >> patchType;
//
//                // Build translation table to convert star patch to foam patch
//                label patchLabel = mapToFoamPatchId[starRegion];
//                if (patchLabel == -1)
//                {
//                    patchLabel = nPatches;
//                    mapToFoamPatchId[starRegion] = patchLabel;
//                    origRegion[patchLabel] = starRegion;
//                    patchTypes_[patchLabel] = patchType;
//
//                    maxId = max(maxId, starRegion);
//
//                    // should actually be case-insensitive
//                    if (patchType == "BAFF")
//                    {
//                        ++nBafflePatches;
//                    }
//                    ++nPatches;
//                }
//
//                ++nPatchFaces[patchLabel];
//            }
//
//            if (nPatches == 0)
//            {
//                Info<< "No boundary faces in file " << inputName << endl;
//            }
//        }
//        else
//        {
//            Info<< "Could not read boundary file " << inputName << endl;
//        }
//    }
//
//    // keep empty patch region in reserve
//    ++nPatches;
//    Info<< "Number of patches = " << nPatches
//        << " (including extra for missing)" << endl;
//
//    // resize
//    origRegion.setSize(nPatches);
//    patchTypes_.setSize(nPatches);
//    patchNames_.setSize(nPatches);
//    nPatchFaces.setSize(nPatches);
//
//    // add our empty patch
//    origRegion[nPatches-1] = 0;
//    nPatchFaces[nPatches-1] = 0;
//    patchTypes_[nPatches-1] = "none";
//
//    // create names
//    // - use 'Label' entry from "constant/boundaryRegion" dictionary
//    forAll(patchTypes_, patchi)
//    {
//        bool fndName = false, fndType = false;
//
//        auto iter = boundaryRegion_.cfind(origRegion[patchi]);
//
//        if (iter.found())
//        {
//            const dictionary& dict = *iter;
//
//            fndType = dict.readIfPresent("BoundaryType", patchTypes_[patchi]);
//            fndName = dict.readIfPresent("Label", patchNames_[patchi]);
//        }
//
//        // Consistent names. Long form and in lowercase
//        if (!fndType)
//        {
//            stringOps::inplaceLower(patchTypes_[patchi]);
//
//            if (patchTypes_[patchi] == "symp")
//            {
//                patchTypes_[patchi] = "symplane";
//            }
//            else if (patchTypes_[patchi] == "cycl")
//            {
//                patchTypes_[patchi] = "cyclic";
//            }
//            else if (patchTypes_[patchi] == "baff")
//            {
//                patchTypes_[patchi] = "baffle";
//            }
//            else if (patchTypes_[patchi] == "moni")
//            {
//                patchTypes_[patchi] = "monitoring";
//            }
//        }
//
//        // Create a name if needed
//        if (!fndName)
//        {
//            patchNames_[patchi] =
//                patchTypes_[patchi] + "_" + name(origRegion[patchi]);
//        }
//    }
//
//    // Enforce name "Default_Boundary_Region"
//    patchNames_[nPatches-1] = defaultBoundaryName;
//
//    // Sort according to ascending region numbers, but leave
//    // "Default_Boundary_Region" as the final patch
//    {
//        labelList sortedIndices
//        (
//            sortedOrder(SubList<label>(origRegion, nPatches-1))
//        );
//
//        labelList oldToNew = identity(nPatches);
//        forAll(sortedIndices, i)
//        {
//            oldToNew[sortedIndices[i]] = i;
//        }
//
//        inplaceReorder(oldToNew, origRegion);
//        inplaceReorder(oldToNew, patchTypes_);
//        inplaceReorder(oldToNew, patchNames_);
//        inplaceReorder(oldToNew, nPatchFaces);
//    }
//
//    // re-sort to have baffles near the end
//    nBafflePatches = 1;
//    if (nBafflePatches)
//    {
//        labelList oldToNew = identity(nPatches);
//        label newIndex = 0;
//        label baffleIndex = (nPatches-1 - nBafflePatches);
//
//        for (label i=0; i < oldToNew.size()-1; ++i)
//        {
//            if (patchTypes_[i] == "baffle")
//            {
//                oldToNew[i] = baffleIndex++;
//            }
//            else
//            {
//                oldToNew[i] = newIndex++;
//            }
//        }
//
//        inplaceReorder(oldToNew, origRegion);
//        inplaceReorder(oldToNew, patchTypes_);
//        inplaceReorder(oldToNew, patchNames_);
//        inplaceReorder(oldToNew, nPatchFaces);
//    }
//
//    mapToFoamPatchId.setSize(maxId+1, -1);
//    forAll(origRegion, patchi)
//    {
//        mapToFoamPatchId[origRegion[patchi]] = patchi;
//    }
//
//    boundaryIds_.setSize(nPatches);
//    forAll(boundaryIds_, patchi)
//    {
//        boundaryIds_[patchi].setSize(nPatchFaces[patchi]);
//        nPatchFaces[patchi] = 0;
//    }
//
//
//    // Pass 2:
//    //
//    if (nPatches > 1 && mapToFoamCellId_.size() > 1)
//    {
//        IFstream is(inputName);
//        readHeader(is, ensightCore::HEADER_BND);
//
//        while (is.read(tok).good() && tok.isLabel())
//        {
//            // Ignore boundary id (not needed)
//
//            is
//                >> starCellId
//                >> cellFaceId
//                >> starRegion
//                >> configNumber
//                >> patchType;
//
//            label patchi = mapToFoamPatchId[starRegion];
//
//            // zero-based indexing
//            cellFaceId--;
//
//            label cellId = -1;
//
//            // convert to FOAM cell number
//            if (starCellId < mapToFoamCellId_.size())
//            {
//                cellId = mapToFoamCellId_[starCellId];
//            }
//
//            if (cellId < 0)
//            {
//                Info
//                    << "Boundaries inconsistent with cell file. "
//                    << "Star cell " << starCellId << " does not exist"
//                    << endl;
//            }
//            else
//            {
//                // restrict lookup to volume cells (no baffles)
//                if (cellId < cellShapes_.size())
//                {
//                    label mapIndex = cellShapes_[cellId].model().index();
//                    if (shapeLookup.found(mapIndex))
//                    {
//                        mapIndex = shapeLookup[mapIndex];
//                        cellFaceId =
//                            ensightCore::starToFoamFaceAddr
//                            [mapIndex][cellFaceId];
//                    }
//                }
//                else
//                {
//                    // we currently use cellId >= nCells to tag baffles,
//                    // we can also use a negative face number
//                    cellFaceId = -1;
//                }
//
//                boundaryIds_[patchi][nPatchFaces[patchi]] =
//                    cellFaceIdentifier(cellId, cellFaceId);
//
//#ifdef DEBUG_BOUNDARY
//                Info<< "bnd " << cellId << " " << cellFaceId << endl;
//#endif
//                // Increment counter of faces in current patch
//                ++nPatchFaces[patchi];
//            }
//        }
//    }
//
//    // retain original information in patchPhysicalTypes_ - overwrite latter
//    patchPhysicalTypes_.setSize(patchTypes_.size());
//
//
//    forAll(boundaryIds_, patchi)
//    {
//        // resize - avoid invalid boundaries
//        if (nPatchFaces[patchi] < boundaryIds_[patchi].size())
//        {
//            boundaryIds_[patchi].setSize(nPatchFaces[patchi]);
//        }
//
//        word origType = patchTypes_[patchi];
//        patchPhysicalTypes_[patchi] = origType;
//
//        if (origType == "symplane")
//        {
//            patchTypes_[patchi] = symmetryPolyPatch::typeName;
//            patchPhysicalTypes_[patchi] = patchTypes_[patchi];
//        }
//        else if (origType == "wall")
//        {
//            patchTypes_[patchi] = wallPolyPatch::typeName;
//            patchPhysicalTypes_[patchi] = patchTypes_[patchi];
//        }
//        else if (origType == "cyclic")
//        {
//            // incorrect. should be cyclicPatch but this
//            // requires info on connected faces.
//            patchTypes_[patchi] = oldCyclicPolyPatch::typeName;
//            patchPhysicalTypes_[patchi] = patchTypes_[patchi];
//        }
//        else if (origType == "baffle")
//        {
//            // incorrect. tag the patch until we get proper support.
//            // set physical type to a canonical "baffle"
//            patchTypes_[patchi] = emptyPolyPatch::typeName;
//            patchPhysicalTypes_[patchi] = "baffle";
//        }
//        else
//        {
//            patchTypes_[patchi] = polyPatch::typeName;
//        }
//
//        Info<< "patch " << patchi
//            << " (region " << origRegion[patchi]
//            << ": " << origType << ") type: '" << patchTypes_[patchi]
//            << "' name: " << patchNames_[patchi] << endl;
//    }
//
//    // cleanup
//    mapToFoamCellId_.clear();
//    cellShapes_.clear();
//}
//
//
////
//// remove unused points
////
//void Foam::fileFormats::ensightMeshReader::cullPoints()
//{
//    label nPoints = points_.size();
//    labelList oldToNew(nPoints, -1);
//
//    // loop through cell faces and note which points are being used
//    forAll(cellFaces_, celli)
//    {
//        const faceList& faces = cellFaces_[celli];
//        forAll(faces, i)
//        {
//            const labelList& labels = faces[i];
//            forAll(labels, j)
//            {
//                oldToNew[labels[j]]++;
//            }
//        }
//    }
//
//    // The new ordering and the count of unused points
//    label pointi = 0;
//    forAll(oldToNew, i)
//    {
//        if (oldToNew[i] >= 0)
//        {
//            oldToNew[i] = pointi++;
//        }
//    }
//
//    // Report unused points
//    if (nPoints > pointi)
//    {
//        Info<< "Unused    points  = " << (nPoints - pointi) << endl;
//        nPoints = pointi;
//
//        // Adjust points and truncate
//        inplaceReorder(oldToNew, points_);
//        points_.setSize(nPoints);
//
//        // Adjust cellFaces - with mesh shapes this might be faster
//        for (faceList& faces : cellFaces_)
//        {
//            for (face& f : faces)
//            {
//                inplaceRenumber(oldToNew, f);
//            }
//        }
//
//        // Adjust baffles
//        for (face& f : baffleFaces_)
//        {
//            inplaceRenumber(oldToNew, f);
//        }
//    }
//}

const Foam::face& Foam::fileFormats::ensightMeshReader::rotateFace
(
    const face& f,
    face& rotatedFace
) const
{
    label fp = findMin(f);

    rotatedFace.setSize(f.size());
    forAll(rotatedFace, i)
    {
        rotatedFace[i] = f[fp];
        fp = f.fcIndex(fp);
    }
    return rotatedFace;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::fileFormats::ensightMeshReader::readGeometry
(
    const scalar scaleFactor
)
{
    ensightReadFile is(geometryFile_);

    string header;
    is.read(header);
    Info<< "Ensight : " << header << endl;
    is.read(header);
    Info<< "Ensight : " << header << endl;

    // Work
    DynamicList<label> verts; 

    label partIndex = -1;
    List<string> partNames;

    label internalPartIndex = -1;

    PtrList<pointField> partPoints;
    PtrList<labelList> partMapToFoamPointIds;

    // Addressing from cell to faces
    // Note: only one internal mesh supported.
    faceListList meshCellFaces;
    labelList meshMapToFoamCellId;

    // Boundary faces
    PtrList<faceList> partFaces;

    string key;
    while (is.good())
    {
        is.readKeyword(key);
        if (key == "node" || key == "element")
        {
            word id;
            is.read(id);
            word op;
            is.read(op);
        }
        else if (key == "part")
        {
            if (partIndex != -1)
            {
                Pout<< decrIndent;
            }

            is.read(partIndex);

            if (partIndex >= partNames.size())
            {
                partNames.setSize(partIndex+1);
            }
            is.read(partNames[partIndex]);

            if (partNames[partIndex] == "internalMesh")
            {
                internalPartIndex = partIndex;
            }
            Pout<< indent
                << "Reading part " << partIndex
                << " name " << partNames[partIndex] << incrIndent << endl;
        }
        else if (key == "coordinates")
        {
            label nPoints;
            is.read(nPoints);

            partPoints.setSize(max(partPoints.size(), partIndex+1));
            partMapToFoamPointIds.setSize(max(partPoints.size(), partIndex+1));
            if (!partPoints.set(partIndex))
            {
                partPoints.set(partIndex, new pointField(nPoints));
                partMapToFoamPointIds.set
                (
                    partIndex,
                    new labelList(nPoints+1, -1)
                );
            }
            auto& points = partPoints[partIndex];
            auto& partMapToFoamPointId = partMapToFoamPointIds[partIndex];

            Pout<< indent << "coordinates " << nPoints << endl;
            for (label pointi = 0; pointi < nPoints; pointi++)
            {
                is.read(points[pointi].x());
                partMapToFoamPointId[pointi+1] = pointi;
            }
            for (label pointi = 0; pointi < nPoints; pointi++)
            {
                is.read(points[pointi].y());
            }
            for (label pointi = 0; pointi < nPoints; pointi++)
            {
                is.read(points[pointi].z());
            }
        }
        else if (key == "hexa8")
        {
            label nShapes;
            is.read(nShapes);

            label celli = meshCellFaces.size();
            meshCellFaces.resize(celli+nShapes);
            meshMapToFoamCellId.setSize(meshCellFaces.size());

            const auto& partMapToFoamPointId = partMapToFoamPointIds[partIndex];

            Pout<< indent<< "hexa8 " << nShapes << endl;
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                verts.clear();
                for (label i = 0; i < 8; i++)
                {
                    label verti;
                    is.read(verti);
                    verts.append(partMapToFoamPointId[verti]);
                }
                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::HEX),
                    verts
                );
                meshCellFaces[celli] = cellVerts.faces();
                meshMapToFoamCellId[celli] = celli;

                celli++;
            }
        }
        //else if (key == "pyramid5")
        //{
        //    label nShapes;
        //    is.read(nShapes);
        //
        //    label celli = partShapes.size();
        //    partShapes.resize(celli+nShapes);
        //
        //    Pout<< "pyramid5:" << nShapes << endl;
        //    for (label shapei = 0; shapei < nShapes; shapei++)
        //    {
        //        verts.clear();
        //        for (label i = 0; i < 8; i++)
        //        {
        //            label verti;
        //            is.read(verti);
        //            verts.append(partMapToFoamPointId[verti]);
        //        }
        //        partShapes[celli++] = cellShape
        //        (
        //            *cellModel::ptr(cellModel::PYR),
        //            verts
        //        );
        //    }
        //}
        else if (key == "nfaced")
        {
            label nShapes;
            is.read(nShapes);

            label celli = meshCellFaces.size();
            meshCellFaces.resize(celli+nShapes);
            meshMapToFoamCellId.setSize(meshCellFaces.size());

            const auto& partMapToFoamPointId = partMapToFoamPointIds[partIndex];

            Pout<< indent<< "nfaces " << nShapes << endl;
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                label nFaces;
                is.read(nFaces);

                faceList& cellFaces = meshCellFaces[celli];
                cellFaces.setSize(nFaces);
                meshMapToFoamCellId[celli] = celli;

                forAll(cellFaces, cellFacei)
                {
                    label nVerts;
                    is.read(nVerts);
                    cellFaces[cellFacei].setSize(nVerts);
                }

                forAll(cellFaces, cellFacei)
                {
                    face& f = cellFaces[cellFacei];
                    forAll(f, fp)
                    {
                        label verti;
                        is.read(verti);
                        f[fp] = partMapToFoamPointId[verti];
                    }
                }

                celli++;
            }
        }
        else if (key == "quad4")
        {
            label nShapes;
            is.read(nShapes);
            Pout<< indent << "quad4 " << nShapes << endl;

            partFaces.setSize(max(partFaces.size(), partIndex+1));
            if (!partFaces.set(partIndex))
            {
                partFaces.set(partIndex, new faceList(0));
            }
            auto& faces = partFaces[partIndex];
            label facei = faces.size();
            faces.setSize(facei+nShapes);

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei];
                f.setSize(4);
                forAll(f, fp)
                {
                    label verti;
                    is.read(verti);
                    f[fp] = verti-1;
                }

                facei++;
            }
        }
        else if (key == "nsided")
        {
            label nShapes;
            is.read(nShapes);
            Pout<< indent << "nsided " << nShapes << endl;

            partFaces.setSize(max(partFaces.size(), partIndex+1));
            if (!partFaces.set(partIndex))
            {
                partFaces.set(partIndex, new faceList(0));
            }
            auto& faces = partFaces[partIndex];
            label facei = faces.size();
            faces.setSize(facei+nShapes);

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                label nVerts;
                is.read(nVerts);
                auto& f = faces[facei];
                f.setSize(nVerts);
                forAll(f, fp)
                {
                    label verti;
                    is.read(verti);
                    f[fp] = verti-1;
                }

                facei++;
            }
        }
    }

    if (partIndex != -1)
    {
        Pout<< decrIndent;
    }


    // Transfer internal points/cellFaces to storage
    points_ = partPoints[internalPartIndex];
    mapToFoamPointId_ = partMapToFoamPointIds[internalPartIndex];

    cellFaces_.transfer(meshCellFaces);
    mapToFoamCellId_.transfer(meshMapToFoamCellId);
    origCellId_ = mapToFoamCellId_;

    // Build map from face to cell and index. Note: use exact match
    // - no circular equivalence
    // but instead pass in ordered faces (lowest numbered vertex first)
    HashTable<cellFaceIdentifier, face, face::symmHasher> vertsToCell;
    {
        label nCellFaces = 0;
        forAll(cellFaces_, celli)
        {
            const auto& cFaces = cellFaces_[celli];
            nCellFaces += cFaces.size();
        }
        vertsToCell.resize(nCellFaces);
    }


    // Insert cell's faces into hashtable
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    face rotatedFace;
    forAll(cellFaces_, celli)
    {
        const auto& cFaces = cellFaces_[celli];
        forAll(cFaces, cFacei)
        {
            const face& f = rotateFace(cFaces[cFacei], rotatedFace);

            const auto fFnd = vertsToCell.find(f);
            if (fFnd)
            {
                // Already inserted. Internal face.
                //Pout<< "For cell:" << celli
                //    << " face:" << cFaces[cFacei]
                //    << " detected other cell/face:" << fFnd() << endl;
                vertsToCell.erase(fFnd);
            }
            else
            {
                //Pout<< "For cell:" << celli
                //    << " face:" << cFaces[cFacei]
                //    << " inserting cell/face:"
                //    << cellFaceIdentifier(celli, cFacei) << endl;
                vertsToCell.insert(f, cellFaceIdentifier(celli, cFacei));
            }
        }
    }


    labelList patchToPart(partNames.size());

    label nPatches = 0;
    forAll(partFaces, parti)
    {
        if (parti != internalPartIndex && partFaces.set(parti))
        {
            patchToPart[nPatches++] = parti;
        }
    }
    patchToPart.setSize(nPatches);

    boundaryIds_.setSize(nPatches);
    patchTypes_.setSize(nPatches, "wall");
    patchNames_.setSize(nPatches);
    forAll(patchNames_, patchi)
    {
        const label parti = patchToPart[patchi];

        Pout<< "Matching part " << parti
            << " name " << partNames[parti]
            << " points " << partPoints[parti].size()
            << " to internal points " << partPoints[internalPartIndex].size()
            << endl;

        patchNames_[patchi] = partNames[parti];

        // Match points to mesh points
        const auto& pp = partPoints[parti];

        labelList partToMesh;
        bool ok = matchPoints
        (
            pp,
            partPoints[internalPartIndex],
            scalarField(pp.size(), ROOTSMALL),
            true,   //false,
            partToMesh
        );

        if (!ok)
        {
            FatalErrorInFunction << "No full match for patch "
                << patchNames_[patchi]  << exit(FatalError);
        }

        const auto& faces = partFaces[parti];

        auto& partCellAndFace = boundaryIds_[patchi];
        partCellAndFace.setSize(faces.size());
        forAll(faces, facei)
        {
            // Rewrite into mesh-point ordering
            const face newF(partToMesh, faces[facei]);
            // Lookup cell and face on cell using the vertices
            const auto cAndF = vertsToCell.find(rotateFace(newF, rotatedFace));
            partCellAndFace[facei] = cAndF();
            vertsToCell.erase(cAndF);
        }
    }
    patchPhysicalTypes_.setSize(nPatches, "wall");
    patchStarts_.setSize(nPatches, 0);
    patchSizes_.setSize(nPatches, 0);

    if (vertsToCell.size())
    {
        // Unused internal or boundary faces
        DebugVar(vertsToCell);
        boundaryIds_.append(List<cellFaceIdentifier>(0));
        {
            auto& cellAndFaces = boundaryIds_.last();
            cellAndFaces.setSize(vertsToCell.size());
            label i = 0;
            for (const auto& e : vertsToCell)
            {
                cellAndFaces[i++] = e;
            }
        }
        patchTypes_.append("empty");
        patchNames_.append("defaultFaces");
        patchPhysicalTypes_.append("empty");
        patchStarts_.append(0);
        patchSizes_.append(0);

        Pout<< "Introducing patch " << patchNames_.last()
            << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::ensightMeshReader::ensightMeshReader
(
    const fileName& prefix,
    const objectRegistry& registry,
    const scalar scaleFactor,
    const bool keepSolids
)
:
    meshReader(prefix, scaleFactor),
    keepSolids_(keepSolids),
    //cellShapes_(0),
    mapToFoamPointId_(0),
    mapToFoamCellId_(0)
{
//    readAux(registry);
}


// ************************************************************************* //
