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

    string key;
    while (is.good())
    {
        do
        {
            is.readKeyword(key);
        }
        while (key.empty() && is.good());

        if (key == "part")
        {
            return false;
        }
        //else if (key == "extents")
        //{
        //    is.read(extents.min().x());
        //    is.read(extents.max().x());
        //    is.read(extents.min().y());
        //    is.read(extents.max().y());
        //    is.read(extents.min().z());
        //    is.read(extents.max().z());
        //    Pout<< indent << "Read extents " << boundBox(min, max) << endl;
        //}
        //else if (key == "node")
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
        //else if (key == "element")
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
        else if (key == "node_ids")
        {
            const label nPoints = points.size();
            // Ignore for now
            for (label i = 0; i < nPoints; i++)
            {
                label index;
                is.read(index);
            }
        }
        else if (key == "coordinates")
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
        else if (key == "tetra4")
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
        else if (key == "pyramid5")
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
                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::PYR),
                    verts
                );
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (key == "penta6")
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
        else if (key == "hexa8")
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
        else if (key == "nfaced")
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
        else if (key == "tria3")
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
        else if (key == "quad4")
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
        else if (key == "nsided")
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
//XXXXXXX


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

    // Work
    label partIndex = -1;
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
    string key;
    while (is.good())
    {
        do
        {
            is.readKeyword(key);
        }
        while (key.empty() && is.good());

        if (key == "extents")
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
        else if (key == "node")
        {
            word id;
            is.read(id);
            word op;
            is.read(op);    // 'none' or 'assign'
            if (op == "given" || op == "ignore")
            {
                read_node_ids = true;
            }
        }
        else if (key == "element")
        {
            word id;
            is.read(id);
            word op;
            is.read(op);    // 'none' or 'assign'
            if (op == "given" || op == "ignore")
            {
                read_elem_ids = true;
            }
        }
        else if (key == "part")
        {
            label partIndex;
            is.read(partIndex);

            bool finished = false;
            do
            {
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
                    << "Reading part " << partIndex
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
                Pout<< decrIndent;
            }
            while (!finished);

            break;
        }
    }



    // Extract the meshPartIndex part
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    // Transfer internal points/cellFaces to storage
    points_ = partPoints[meshPartIndex_];
    //mapToFoamPointId_ = partMapToFoamPointIds[meshPartIndex_];

    cellFaces_.transfer(partCells[meshPartIndex_]);
    //mapToFoamCellId_.transfer(meshMapToFoamCellId);
    //origCellId_ = mapToFoamCellId_;

    // Build map from face to cell and index. Note: use exact match
    // - no circular equivalence
    // - but instead pass in ordered faces (lowest numbered vertex first)
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
                vertsToCell.erase(fFnd);
            }
            else
            {
                vertsToCell.insert(f, cellFaceIdentifier(celli, cFacei));
            }
        }
    }


    DynamicList<label> patchParts;
    if (useAllPartsForPatching_)
    {
        forAll(partFaces, parti)
        {
            if (partFaces[parti].size())
            {
                patchParts.append(parti);
            }
        }
        Pout<< "Using parts " << UIndirectList<string>(partNames, patchParts)
            << " for patching" << endl;
    }
    else
    {
        patchParts.append(meshPartIndex_);
        Pout<< "Using part " << partNames[meshPartIndex_]
            << " for patching" << endl;
    }


    labelList patchToPart(partNames.size());

    label nPatches = 0;

    for (const label parti : patchParts)
    {
        patchToPart[nPatches++] = parti;
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
            << " to internal points " << partPoints[meshPartIndex_].size()
            << endl;

        patchNames_[patchi] = partNames[parti];

        // Match points to mesh points
        const auto& pp = partPoints[parti];

        labelList partToMesh;
        if (parti != meshPartIndex_)
        {
            bool ok = matchPoints
            (
                pp,
                partPoints[meshPartIndex_],
                scalarField(pp.size(), ROOTSMALL),
                true,   //false,
                partToMesh
            );

            if (!ok)
            {
                FatalErrorInFunction << "No full match for patch "
                    << patchNames_[patchi]  << exit(FatalError);
            }
        }
        else
        {
            partToMesh = identity(pp.size());
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
    const label meshPartIndex,
    const bool useAllPartsForPatching,
    const objectRegistry& registry,
    const scalar scaleFactor
    //const bool keepSolids
)
:
    meshReader(geomFile, scaleFactor),
    meshPartIndex_(meshPartIndex),
    useAllPartsForPatching_(useAllPartsForPatching)
    //mapToFoamPointId_(0),
    //mapToFoamCellId_(0)
{}


// ************************************************************************* //
