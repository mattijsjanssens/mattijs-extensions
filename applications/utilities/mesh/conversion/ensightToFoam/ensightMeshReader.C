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
        do
        {
            is.readKeyword(key);
        }
        while (key.empty() && is.good());

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
                << " name " << partNames[partIndex]
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg()
                << incrIndent << endl;
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

            Pout<< indent << "coordinates " << nPoints
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;
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
        else if (key == "tetra4")
        {
            label nShapes;
            is.read(nShapes);

            const label celli = meshCellFaces.size();
            meshCellFaces.resize(celli+nShapes);
            meshMapToFoamCellId.setSize(meshCellFaces.size());

            const auto& partMapToFoamPointId = partMapToFoamPointIds[partIndex];

            Pout<< indent<< "tetra4 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                verts.clear();
                for (label i = 0; i < 4; i++)
                {
                    label verti;
                    is.read(verti);
                    verts.append(partMapToFoamPointId[verti]);
                }
                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::TET),
                    verts
                );
                meshCellFaces[celli+shapei] = cellVerts.faces();
                meshMapToFoamCellId[celli+shapei] = celli+shapei;
            }
        }
        else if (key == "pyramid5")
        {
            label nShapes;
            is.read(nShapes);

            const label celli = meshCellFaces.size();
            meshCellFaces.resize(celli+nShapes);
            meshMapToFoamCellId.setSize(meshCellFaces.size());

            const auto& partMapToFoamPointId = partMapToFoamPointIds[partIndex];

            Pout<< indent<< "pyramid5 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                verts.clear();
                for (label i = 0; i < 5; i++)
                {
                    label verti;
                    is.read(verti);
                    verts.append(partMapToFoamPointId[verti]);
                }
                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::PYR),
                    verts
                );
                meshCellFaces[celli+shapei] = cellVerts.faces();
                meshMapToFoamCellId[celli+shapei] = celli+shapei;
            }
        }
        else if (key == "penta6")
        {
            label nShapes;
            is.read(nShapes);

            const label celli = meshCellFaces.size();
            meshCellFaces.resize(celli+nShapes);
            meshMapToFoamCellId.setSize(meshCellFaces.size());

            const auto& partMapToFoamPointId = partMapToFoamPointIds[partIndex];

            Pout<< indent<< "penta6 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                verts.clear();
                for (label i = 0; i < 6; i++)
                {
                    label verti;
                    is.read(verti);
                    verts.append(partMapToFoamPointId[verti]);
                }
                cellShape cellVerts
                (
                    *cellModel::ptr(cellModel::PRISM),
                    verts
                );
                meshCellFaces[celli+shapei] = cellVerts.faces();
                meshMapToFoamCellId[celli+shapei] = celli+shapei;
            }
        }
        else if (key == "hexa8")
        {
            label nShapes;
            is.read(nShapes);

            const label celli = meshCellFaces.size();
            meshCellFaces.resize(celli+nShapes);
            meshMapToFoamCellId.setSize(meshCellFaces.size());

            const auto& partMapToFoamPointId = partMapToFoamPointIds[partIndex];

            Pout<< indent<< "hexa8 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;
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
                meshCellFaces[celli+shapei] = cellVerts.faces();
                meshMapToFoamCellId[celli+shapei] = celli+shapei;
            }
        }
        else if (key == "nfaced")
        {
            label nShapes;
            is.read(nShapes);

            const label celli = meshCellFaces.size();
            meshCellFaces.resize(celli+nShapes);
            meshMapToFoamCellId.setSize(meshCellFaces.size());

            const auto& partMapToFoamPointId = partMapToFoamPointIds[partIndex];

            Pout<< indent<< "nfaced " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                label nFaces;
                is.read(nFaces);
                faceList& cellFaces = meshCellFaces[celli+shapei];
                cellFaces.setSize(nFaces);
                meshMapToFoamCellId[celli+shapei] = celli+shapei;
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                faceList& cellFaces = meshCellFaces[celli+shapei];
                forAll(cellFaces, cellFacei)
                {
                    label nVerts;
                    is.read(nVerts);
                    cellFaces[cellFacei].setSize(nVerts);
                }
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                faceList& cellFaces = meshCellFaces[celli+shapei];
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
            }
        }
        else if (key == "tria3")
        {
            label nShapes;
            is.read(nShapes);
            Pout<< indent << "tria3 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            partFaces.setSize(max(partFaces.size(), partIndex+1));
            if (!partFaces.set(partIndex))
            {
                partFaces.set(partIndex, new faceList(0));
            }
            auto& faces = partFaces[partIndex];
            const label facei = faces.size();
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

            partFaces.setSize(max(partFaces.size(), partIndex+1));
            if (!partFaces.set(partIndex))
            {
                partFaces.set(partIndex, new faceList(0));
            }
            auto& faces = partFaces[partIndex];
            const label facei = faces.size();
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

            partFaces.setSize(max(partFaces.size(), partIndex+1));
            if (!partFaces.set(partIndex))
            {
                partFaces.set(partIndex, new faceList(0));
            }
            auto& faces = partFaces[partIndex];
            const label facei = faces.size();
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
    const fileName& prefix,
    const objectRegistry& registry,
    const scalar scaleFactor,
    const bool keepSolids
)
:
    meshReader(prefix, scaleFactor),
    mapToFoamPointId_(0),
    mapToFoamCellId_(0)
{}


// ************************************************************************* //
