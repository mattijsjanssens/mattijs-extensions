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

#include "uFieldReconstructor.H"
#include "unallocatedFvMesh.H"
#include "unallocatedFvMeshTools.H"
#include "unallocatedFvFieldReconstructor.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(uFieldReconstructor, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::uFieldReconstructor::readProcDatabases
(
    const IOobject& io,
    const label nProcs
) const
{
    procDatabases_.setSize(nProcs);
    forAll(procDatabases_, proci)
    {
        procDatabases_.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                io.rootPath(),
                io.caseName()/fileName(word("processor") + Foam::name(proci))
            )
        );
    }
}


void Foam::uFieldReconstructor::readProcMeshes
(
    const Time& baseTime,
    const fileName& instance
) const
{
    cellProcAddressing_.setSize(procDatabases_.size());
    faceProcAddressing_.setSize(procDatabases_.size());
    boundaryProcAddressing_.setSize(procDatabases_.size());

    forAll(cellProcAddressing_, proci)
    {
        cellProcAddressing_.set
        (
            proci,
            new labelIOList
            (
                IOobject
                (
                    "cellProcAddressing",
                    instance,
                    fvMesh::meshSubDir,
                    procDatabases_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
    forAll(faceProcAddressing_, proci)
    {
        faceProcAddressing_.set
        (
            proci,
            new labelIOList
            (
                IOobject
                (
                    "faceProcAddressing",
                    instance,
                    fvMesh::meshSubDir,
                    procDatabases_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
    forAll(boundaryProcAddressing_, proci)
    {
        boundaryProcAddressing_.set
        (
            proci,
            new labelIOList
            (
                IOobject
                (
                    "boundaryProcAddressing",
                    instance,
                    fvMesh::meshSubDir,
                    procDatabases_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }

    // Read the (unallocated) processor meshes
    procMeshes_.setSize(procDatabases_.size());

    forAll(procMeshes_, proci)
    {
        procMeshes_.set
        (
            proci,
            unallocatedFvMeshTools::newMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    instance,
                    procDatabases_[proci],
                    IOobject::MUST_READ
                ),
                cellProcAddressing_[proci].size()
            )
        );
        if (debug)
        {
            Pout<< "Read mesh for processor " << proci << ':' << nl
                << incrIndent << procMeshes_[proci].info() << decrIndent
                << endl;
        }
    }


    // Get mesh as unallocated
    baseMeshPtr_ = unallocatedFvMeshTools::newMesh
    (
        IOobject
        (
            fvMesh::defaultRegion,      // name of mesh
            instance,                   // start search
            baseTime,   //baseDatabasePtr_(),
            IOobject::MUST_READ
        )
    );
    if (debug)
    {
        Pout<< "Read baseMesh:" << nl
            << incrIndent << baseMeshPtr_().info() << decrIndent << endl;
    }

    reconstructorPtr_.reset
    (
        new unallocatedFvFieldReconstructor
        (
            baseMeshPtr_(),
            procMeshes_,
            faceProcAddressing_,
            cellProcAddressing_,
            boundaryProcAddressing_
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uFieldReconstructor::uFieldReconstructor(const polyMesh& mesh)
:
    MeshObject<polyMesh, TopologicalMeshObject, uFieldReconstructor>(mesh)
{
   // Read the processor databases
    label nProcs = fileHandler().nProcs(mesh.time().path(), word::null);
    if (debug)
    {
        Pout<< "uFieldReconstructor detected nProcs:" << nProcs
            << " from:" << mesh.time().path() << endl;
    }
    readProcDatabases(mesh.thisDb(), nProcs);
    readProcMeshes(mesh.time(), mesh.time().constant());
}


// ************************************************************************* //
