/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "pointMesh.H"
#include "globalMeshData.H"
#include "pointMeshMapper.H"
#include "pointFields.H"
#include "MapGeometricFields.H"
#include "MapPointField.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointMesh, 0);
    word pointMesh::meshSubDir = "pointMesh";
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMesh::mapFields(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Pout<< "void pointMesh::mapFields(const mapPolyMesh&): "
            << "Mapping all registered pointFields."
            << endl;
    }
    // Create a mapper
    const pointMeshMapper m(*this, mpm);

    MapGeometricFields<scalar, pointPatchField, pointMeshMapper, pointMesh>(m);
    MapGeometricFields<vector, pointPatchField, pointMeshMapper, pointMesh>(m);
    MapGeometricFields
    <
        sphericalTensor,
        pointPatchField,
        pointMeshMapper,
        pointMesh
    >(m);
    MapGeometricFields<symmTensor, pointPatchField, pointMeshMapper, pointMesh>
    (m);
    MapGeometricFields<tensor, pointPatchField, pointMeshMapper, pointMesh>(m);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMesh::pointMesh(const polyMesh& pMesh)
:
    MeshObject<polyMesh, Foam::UpdateableMeshObject, pointMesh>(pMesh),
    GeoMesh<polyMesh>(pMesh),
    meshEdges_
    (
        IOobject
        (
            "meshEdges",
            pMesh.facesInstance(),
            meshSubDir,
            thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        edgeList(0)
    ),
    meshPoints_
    (
        IOobject
        (
            "meshPoints",
            pMesh.facesInstance(),
            meshSubDir,
            thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelList(0)
    ),
    boundary_(*this, pMesh.boundaryMesh())
{
    if (debug)
    {
        Pout<< "pointMesh::pointMesh(const polyMesh&): "
            << "Constructing from polyMesh " << pMesh.name()
            << endl;
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


Foam::pointMesh::pointMesh(const IOobject& io, const polyMesh& pMesh)
:
    MeshObject<polyMesh, Foam::UpdateableMeshObject, pointMesh>(pMesh),
    GeoMesh<polyMesh>(pMesh),
    meshEdges_
    (
        IOobject
        (
            "meshEdges",
            io.instance(),
            meshSubDir,
            thisDb(),
            io.readOpt(),
            io.writeOpt(),
            false
        ),
        edgeList(0)
    ),
    meshPoints_
    (
        IOobject
        (
            "meshPoints",
            io.instance(),
            meshSubDir,
            thisDb(),
            io.readOpt(),
            io.writeOpt(),
            false
        ),
        labelList(0)
    ),
    boundary_(io, *this, pMesh.boundaryMesh())
{
    if (debug)
    {
        Pout<< "pointMesh::pointMesh(const polyMesh&): "
            << "Constructing from IO " << io.objectPath()
            << endl;
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMesh::~pointMesh()
{
    if (debug)
    {
        Pout<< "~pointMesh::pointMesh()"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::pointMesh::meshDir() const
{
    return mesh().dbDir()/meshSubDir;
}


bool Foam::pointMesh::movePoints()
{
    if (debug)
    {
        Pout<< "pointMesh::movePoints(const pointField&): "
            << "Moving points." << endl;
    }

    boundary_.movePoints(GeoMesh<polyMesh>::mesh_.points());

    return true;
}


void Foam::pointMesh::updateMesh(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Pout<< "pointMesh::updateMesh(const mapPolyMesh&): "
            << "Updating for topology changes." << endl;
        Pout<< endl;
    }
    boundary_.updateMesh();

    forAll(meshEdges_, edgei)
    {
        edge& e = meshEdges_[edgei];
        e[0] = mpm.reversePointMap()[e[0]];
        e[1] = mpm.reversePointMap()[e[1]];
    }

    forAll(meshPoints_, pointi)
    {
        meshPoints_[pointi] = mpm.reversePointMap()[meshPoints_[pointi]];
    }

    // Map all registered point fields
    mapFields(mpm);
}


// ************************************************************************* //
