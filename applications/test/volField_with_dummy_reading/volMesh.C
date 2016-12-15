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

\*---------------------------------------------------------------------------*/

#include "volMesh.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMesh> Foam::volMesh::New
(
    const IOobject& io,
    const fvMesh& mesh,
    const bool validBoundary
)
{
    pointField newPoints(0);
    faceList newFaces(0);
    cellList newCells(0);

    autoPtr<fvMesh> dummyMeshPtr
    (
        new fvMesh
        (
            io,
            xferMove(newPoints),
            xferMove(newFaces),
            xferMove(newCells),
            validBoundary               // parallel sync?
        )
    );
    fvMesh& dummyMesh = dummyMeshPtr();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    List<polyPatch*> newBoundary(bm.size());
    forAll(newBoundary, patchi)
    {
        newBoundary[patchi] = bm[patchi].clone
        (
            dummyMesh.boundaryMesh(),
            patchi,
            labelList(0),
            0
        ).ptr();
    }

    dummyMesh.addFvPatches(newBoundary, validBoundary);

    return dummyMeshPtr;
}


Foam::autoPtr<Foam::fvMesh> Foam::volMesh::New
(
    const IOobject& io,
    Istream& is,
    const bool validBoundary
)
{
    pointField newPoints(is);
    faceList newFaces(is);
    cellList newCells(is);
    PtrList<entry> patchEntries(is);

    autoPtr<fvMesh> meshPtr
    (
        new fvMesh
        (
            io,
            xferMove(newPoints),
            xferMove(newFaces),
            xferMove(newCells),
            validBoundary               // parallel sync?
        )
    );



    fvMesh& mesh = meshPtr();

    List<polyPatch*> patches(patchEntries.size());

    label nPatches = 0;
    forAll(patchEntries, patchi)
    {
        const entry& e = patchEntries[patchi];
        const word type(e.dict().lookup("type"));
        const word& name = e.keyword();

        dictionary patchDict(e.dict());
        patchDict.set("nFaces", 0);
        patchDict.set("startFace", 0);

        patches[patchi] = polyPatch::New
        (
            name,
            patchDict,
            nPatches++,
            mesh.boundaryMesh()
        ).ptr();
    }

    patches.setSize(nPatches);
    mesh.addFvPatches(patches, validBoundary);

    return meshPtr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
