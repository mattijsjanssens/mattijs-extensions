/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "unallocatedFvMeshObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(unallocatedFvMeshObject, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::unallocatedFvMeshObject::patchNames
(
    const fvMesh& mesh
)
{
    const fvBoundaryMesh& bm = mesh.boundary();
    wordList names(bm.size());
    forAll(bm, i)
    {
        names[i] = bm[i].name();
    }
    return names;
}


Foam::wordList Foam::unallocatedFvMeshObject::patchTypes
(
    const fvMesh& mesh
)
{
    const fvBoundaryMesh& bm = mesh.boundary();
    wordList types(bm.size());
    forAll(bm, i)
    {
        types[i] = bm[i].type();
    }
    return types;
}


Foam::labelList Foam::unallocatedFvMeshObject::patchSizes
(
    const fvMesh& mesh
)
{
    const fvBoundaryMesh& bm = mesh.boundary();
    labelList sizes(bm.size());
    forAll(bm, i)
    {
        sizes[i] = bm[i].size();
    }
    return sizes;
}


Foam::labelList Foam::unallocatedFvMeshObject::patchStarts
(
    const fvMesh& mesh
)
{
    const fvBoundaryMesh& bm = mesh.boundary();
    labelList starts(bm.size());
    forAll(bm, i)
    {
        starts[i] = bm[i].start();
    }
    return starts;
}


Foam::List<Foam::wordList> Foam::unallocatedFvMeshObject::patchGroups
(
    const fvMesh& mesh
)
{
    const fvBoundaryMesh& bm = mesh.boundary();
    List<wordList> groups(bm.size());
    forAll(bm, i)
    {
        groups[i] = bm[i].patch().inGroups();
    }
    return groups;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unallocatedFvMeshObject::unallocatedFvMeshObject(const fvMesh& mesh)
:
    MeshObject<fvMesh, TopologicalMeshObject, unallocatedFvMeshObject>(mesh),
    unallocatedFvMesh
    (
        mesh.name(),
        mesh.thisDb(),
        mesh.nCells(),
        patchNames(mesh),
        patchTypes(mesh),
        patchSizes(mesh),
        patchStarts(mesh),
        patchGroups(mesh),
        mesh.globalData()
    )
{}


// ************************************************************************* //
