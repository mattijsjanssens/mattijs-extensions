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
#include "cellModel.H"
#include "ensightReadFile.H"
#include "matchPoints.H"
#include "mergePoints.H"
#include "ListListOps.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{
    defineTypeNameAndDebug(ensightMeshReader, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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


bool Foam::fileFormats::ensightMeshReader::readGoldPart
(
    ensightReadFile& is,
    const bool read_node_ids,
    const bool read_elem_ids,

    pointField& points,
    labelList& pointToNodeIds,

    // 3D-elems : cells (cell-to-faces)
    faceListList& cells,
    labelList& cellToElemIds,

    // 2D-elems : faces
    faceList& faces,
    labelList& faceToElemIDs
) const
{
    //- Read a single part. Return true if end-of-file reached. Return false
    //  if 'part' read


    // Work
    DynamicList<label> verts; 

    string line;
    while (is.good())
    {
        do
        {
            is.readKeyword(line);
        }
        while (line.empty() && is.good());
        const auto split = stringOps::splitSpace(line);

        if (split[0] == "part")
        {
            return false;
        }
        //else if (split[0] == "extents")
        //{
        //    is.read(extents.min().x());
        //    is.read(extents.max().x());
        //    is.read(extents.min().y());
        //    is.read(extents.max().y());
        //    is.read(extents.min().z());
        //    is.read(extents.max().z());
        //    Pout<< indent << "Read extents " << boundBox(min, max) << endl;
        //}
        //else if (split[0] == "node")
        //{
        //    word id;
        //    is.read(id);
        //    word op;
        //    is.read(op);    // 'none' or 'assign' or 'given'
        //    if (op == "given" || op == "ignore")
        //    {
        //        read_node_ids = true;
        //    }
        //}
        //else if (split[0] == "element")
        //{
        //    word id;
        //    is.read(id);
        //    word op;
        //    is.read(op);    // 'none' or 'assign' or 'given'
        //    if (op == "given" || op == "ignore")
        //    {
        //        read_elem_ids = true;
        //    }
        //}
        else if (split[0] == "node_ids")
        {
            const label nPoints = points.size();
            // Ignore for now
            for (label i = 0; i < nPoints; i++)
            {
                label index;
                is.read(index);
            }
        }
        else if (split[0] == "coordinates")
        {
            label nPoints;
            is.read(nPoints);

            Pout<< indent << "coordinates " << nPoints
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            if (read_node_ids)
            {
                // Ensight Gold: read node ids separate
                pointToNodeIds.setSize(nPoints, -1);
                for (label pointi = 0; pointi < nPoints; pointi++)
                {
                    is.read(pointToNodeIds[pointi]);
                }
            }

            points.setSize(nPoints);

            for (label pointi = 0; pointi < nPoints; pointi++)
            {
                is.read(points[pointi].x());
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
        else if (split[0] == "tetra4")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "tetra4 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            if (read_elem_ids)
            {
                cellToElemIds.resize(celli+nShapes);
                for (label shapei = 0; shapei < nShapes; shapei++)
                {
                    is.read(cellToElemIds[celli+shapei]);
                }
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                verts.clear();
                for (label i = 0; i < 4; i++)
                {
                    label verti;
                    is.read(verti);
                    verts.append(verti-1);
                }
                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::TET),
                    verts
                );
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (split[0] == "pyramid5")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "pyramid5 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            if (read_elem_ids)
            {
                cellToElemIds.resize(celli+nShapes);
                for (label shapei = 0; shapei < nShapes; shapei++)
                {
                    is.read(cellToElemIds[celli+shapei]);
                }
                DebugVar(cellToElemIds);
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                verts.clear();
                for (label i = 0; i < 5; i++)
                {
                    label verti;
                    is.read(verti);
                    verts.append(verti-1);
                }

                DebugVar(verts);

                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::PYR),
                    verts
                );
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (split[0] == "penta6")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "penta6 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            if (read_elem_ids)
            {
                cellToElemIds.resize(celli+nShapes);
                for (label shapei = 0; shapei < nShapes; shapei++)
                {
                    is.read(cellToElemIds[celli+shapei]);
                }
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                verts.clear();
                for (label i = 0; i < 6; i++)
                {
                    label verti;
                    is.read(verti);
                    verts.append(verti-1);
                }
                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::PRISM),
                    verts
                );
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (split[0] == "hexa8")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "hexa8 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            if (read_elem_ids)
            {
                cellToElemIds.resize(celli+nShapes);
                for (label shapei = 0; shapei < nShapes; shapei++)
                {
                    is.read(cellToElemIds[celli+shapei]);
                }
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                verts.clear();
                for (label i = 0; i < 8; i++)
                {
                    label verti;
                    is.read(verti);
                    verts.append(verti-1);
                }
                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::HEX),
                    verts
                );
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (split[0] == "nfaced")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "nfaced " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            if (read_elem_ids)
            {
                cellToElemIds.resize(celli+nShapes);
                for (label shapei = 0; shapei < nShapes; shapei++)
                {
                    is.read(cellToElemIds[celli+shapei]);
                }
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                label nFaces;
                is.read(nFaces);
                faceList& cellFaces = cells[celli+shapei];
                cellFaces.setSize(nFaces);
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                faceList& cellFaces = cells[celli+shapei];
                forAll(cellFaces, cellFacei)
                {
                    label nVerts;
                    is.read(nVerts);
                    cellFaces[cellFacei].setSize(nVerts);
                }
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                faceList& cellFaces = cells[celli+shapei];
                forAll(cellFaces, cellFacei)
                {
                    face& f = cellFaces[cellFacei];
                    forAll(f, fp)
                    {
                        label verti;
                        is.read(verti);
                        f[fp] = verti-1;
                    }
                }
            }
        }
        else if (split[0] == "tria3")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent << "tria3 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label facei = faces.size();

            if (read_elem_ids)
            {
                faceToElemIDs.resize(facei+nShapes);
                for (label shapei = 0; shapei < nShapes; shapei++)
                {
                    is.read(faceToElemIDs[facei+shapei]);
                }
            }

            faces.setSize(facei+nShapes);

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei+shapei];
                f.setSize(3);
                forAll(f, fp)
                {
                    label verti;
                    is.read(verti);
                    f[fp] = verti-1;
                }
            }
        }
        else if (split[0] == "quad4")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent << "quad4 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label facei = faces.size();

            if (read_elem_ids)
            {
                faceToElemIDs.resize(facei+nShapes);
                for (label shapei = 0; shapei < nShapes; shapei++)
                {
                    is.read(faceToElemIDs[facei+shapei]);
                }
            }

            faces.setSize(facei+nShapes);

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei+shapei];
                f.setSize(4);
                forAll(f, fp)
                {
                    label verti;
                    is.read(verti);
                    f[fp] = verti-1;
                }
            }
        }
        else if (split[0] == "nsided")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent << "nsided " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label facei = faces.size();

            if (read_elem_ids)
            {
                faceToElemIDs.resize(facei+nShapes);
                for (label shapei = 0; shapei < nShapes; shapei++)
                {
                    is.read(faceToElemIDs[facei+shapei]);
                }
            }

            faces.setSize(facei+nShapes);

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei+shapei];
                label nVerts;
                is.read(nVerts);
                f.setSize(nVerts);
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei+shapei];
                forAll(f, fp)
                {
                    label verti;
                    is.read(verti);
                    f[fp] = verti-1;
                }
            }
        }
    }

    // EOF
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::fileFormats::ensightMeshReader::readGeometry
(
    const scalar scaleFactor
)
{
    ensightReadFile is(geometryFile_);

    // Skip 'binary' tag
    is.readBinaryHeader();

    string header;
    is.read(header);
    Info<< "Ensight : " << header << endl;
    is.read(header);
    Info<< "Ensight : " << header << endl;


    bool read_node_ids = false;
    bool read_elem_ids = false;

    // Storage for all the parts
    List<string> partNames;
    PtrList<pointField> partPoints;
    PtrList<labelList> partNodeIDs;

    // Cells (cell-to-faces)
    // Note: only one internal mesh supported.
    PtrList<faceListList> partCells;
    // Element ID for cells
    PtrList<labelList> partCellIDs;
    // Boundary faces
    PtrList<faceList> partFaces;
    // Element IDs for faces
    PtrList<labelList> partFaceIDs;


    // Parse all
    string line;
    while (is.good())
    {
        do
        {
            is.readKeyword(line);
        }
        while (line.empty() && is.good());
        const auto split = stringOps::splitSpace(line);

        if (split[0] == "extents")
        {
            point min;
            point max;
            is.read(min.x());
            is.read(max.x());
            is.read(min.y());
            is.read(max.y());
            is.read(min.z());
            is.read(max.z());
            Pout<< indent
                << "Read extents " << boundBox(min, max)
                << endl;
        }
        else if (split[0] == "node")
        {
            word id(split[1]);
            word op(split[2]);
            if (op == "given" || op == "ignore")
            {
                Pout<< indent << "Reading node ids" << endl;
                read_node_ids = true;
            }
        }
        else if (split[0] == "element")
        {
            word id(split[1]);
            word op(split[2]);
            if (op == "given" || op == "ignore")
            {
                Pout<< indent << "Reading element ids" << endl;
                read_elem_ids = true;
            }
        }
        else if (split[0] == "part")
        {
            bool finished = false;
            do
            {
                label partIndex;
                is.read(partIndex);
                partIndex -= 1;

                if (partIndex >= partNames.size())
                {
                    partNames.setSize(partIndex+1);
                    partPoints.setSize(partIndex+1);
                    partPoints.set(partIndex, new pointField(0));
                    partNodeIDs.setSize(partIndex+1);
                    partNodeIDs.set(partIndex, new labelList(0));
                    partCells.setSize(partIndex+1);
                    partCells.set(partIndex, new faceListList(0));
                    partCellIDs.setSize(partIndex+1);
                    partCellIDs.set(partIndex, new labelList(0));
                    partFaces.setSize(partIndex+1);
                    partFaces.set(partIndex, new faceList(0));
                    partFaceIDs.setSize(partIndex+1);
                    partFaceIDs.set(partIndex, new labelList(0));
                }

                is.read(partNames[partIndex]);

                Pout<< indent
                    << "Reading part " << partIndex+1
                    << " name " << partNames[partIndex]
                    << " starting at line " << is.lineNumber()
                    << " position " << is.stdStream().tellg() << endl;

                Pout<< incrIndent;
                finished = readGoldPart
                (
                    is,
                    read_node_ids,
                    read_elem_ids,

                    partPoints[partIndex],
                    partNodeIDs[partIndex],

                    // Cells (cell-to-faces)
                    partCells[partIndex],
                    partCellIDs[partIndex],

                    // Faces
                    partFaces[partIndex],
                    partFaceIDs[partIndex]
                );

                Pout<< indent
                    << "For part " << partIndex+1
                    << " read cells " << partCells[partIndex].size()
                    << " faces " << partFaces[partIndex].size()
                    << endl;

                Pout<< decrIndent;
            }
            while (!finished);

            break;
        }
    }


    // Merge all coordinates
    labelListList pointToMerged(partPoints.size());
    {
        label nPoints = 0;
        forAll(partPoints, parti)
        {
            nPoints += partPoints[parti].size();
        }

        points_.setSize(nPoints);
        nPoints = 0;
        forAll(partPoints, parti)
        {
            const auto& pts = partPoints[parti];
            SubList<point>(points_, pts.size(), nPoints) = pts;
            auto& map = pointToMerged[parti];
            map = nPoints + identity(pts.size());
            nPoints += pts.size();
        }

        if (mergeTol_ > 0)
        {
            labelList sharedToMerged;
            const label nMerged = inplaceMergePoints
            (
                points_,
                mergeTol_,
                false,
                sharedToMerged
            );
            Pout<< "Merged " << nMerged << " points out of " << nPoints
                << endl;

            forAll(partPoints, parti)
            {
                inplaceRenumber(sharedToMerged, pointToMerged[parti]);
            }
        }
    }

    // Merge all cells
    labelList cellOffsets(partCells.size()+1);
    cellOffsets[0] = 0;
    {
        label nCells = 0;
        forAll(partCells, parti)
        {
            nCells += partCells[parti].size();
            cellOffsets[parti+1] = nCells;
        }
        DebugVar(cellOffsets);


        cellFaces_.setSize(nCells);
        cellTableId_.setSize(nCells);
        forAll(partCells, parti)
        {
            // Copy cells into position
            const auto& cells = partCells[parti];

            SubList<faceList> cellSlice
            (
                cellFaces_,
                cells.size(),
                cellOffsets[parti]
            );
            cellSlice = cells;

            Pout<< "part:" << parti
                << " slice:" << cellSlice.size()
                << " start:" << cellOffsets[parti]
                << endl;

            SubList<label> cellIDSlice
            (
                cellTableId_,
                cells.size(),
                cellOffsets[parti]
            );
            cellIDSlice = parti;


            // Renumber faces
            const auto& pointMap = pointToMerged[parti];
            Pout<< "For part:" << parti
                << " using pointMap:" << flatOutput(pointMap) << endl;

            label celli = cellOffsets[parti];
            for (auto& cellFaces : cellSlice)
            {
                //auto& cellFaces = cellFaces_[celli];

                for (auto& f : cellFaces)
                {
                    inplaceRenumber(pointMap, f);
                }
                Pout<< "cell:" << celli
                    << " faces:" << flatOutput(cellFaces) << endl;
                celli++;
            }
        }
    }


    // Add cells to separate zones
    forAll(partCells, parti)
    {
        cellTable_.setMaterial(parti, "fluid");
        cellTable_.setName(parti, "part" + Foam::name(parti));
    }


    // Build map from face to cell and index. Note: use exact match
    // - no circular equivalence
    // - but instead pass in ordered faces (lowest numbered vertex first)
    HashTable<cellFaceIdentifier, face, face::symmHasher> vertsToCell
    (
        2*cellOffsets.last()
    );

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
                vertsToCell.erase(fFnd);
            }
            else
            {
                vertsToCell.insert(f, cellFaceIdentifier(celli, cFacei));
            }
        }
    }


    labelList patchToPart(partNames.size());
    label nPatches = 0;
    forAll(partFaces, parti)
    {
        if (partFaces[parti].size())
        {
            Pout<< "Using part " << parti
                << " name " << partNames[parti]
                << " as patch " << nPatches
                << endl;

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
            << " to merged points " << points_.size()
            << endl;

        patchNames_[patchi] = partNames[parti];

        const auto& pointMap = pointToMerged[parti];
        const auto& faces = partFaces[parti];

        auto& partCellAndFace = boundaryIds_[patchi];
        partCellAndFace.setSize(faces.size());
        forAll(faces, facei)
        {
            // Rewrite into mesh-point ordering
            const face newF(pointMap, faces[facei]);
            // Lookup cell and face on cell using the vertices
            const auto cAndF = vertsToCell.find(rotateFace(newF, rotatedFace));

            if (!cAndF)
            {
                FatalErrorInFunction << "Trying to find face " << facei
                    << " verts:" << newF
                    << " rotatedFace:" << rotatedFace
                    << " in table:" << vertsToCell
                    << exit(FatalError);
            }
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

        Pout<< "Introducing default patch " << patchNames_.size()-1
            << " name " << patchNames_.last() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::ensightMeshReader::ensightMeshReader
(
    const fileName& geomFile,
    const objectRegistry& registry,
    const scalar mergeTol,
    const scalar scaleFactor
)
:
    meshReader(geomFile, scaleFactor),
    mergeTol_(mergeTol)
{}


// ************************************************************************* //
