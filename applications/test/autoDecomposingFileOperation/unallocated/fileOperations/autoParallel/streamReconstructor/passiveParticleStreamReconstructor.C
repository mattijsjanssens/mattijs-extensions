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

#include "passiveParticleStreamReconstructor.H"
#include "uFieldReconstructor.H"
#include "addToRunTimeSelectionTable.H"
#include "unallocatedIOPosition.H"

namespace Foam
{
    defineTypeName(passiveParticleStreamReconstructor);
    addToRunTimeSelectionTable
    (
        streamReconstructor,
        passiveParticleStreamReconstructor,
        cloudName
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::passiveParticleStreamReconstructor::reconstruct
(
    const IOobject& io,
    const bool,
    Ostream& os
) const
{
    Pout<< "*** passiveParticleStreamReconstructor : "
        << io.objectPath() << endl;

    const uFieldReconstructor& reconstructor =
        uFieldReconstructor::New(io.db().parent().parent());

    const PtrList<unallocatedFvMesh>& procMeshes = reconstructor.procMeshes();

    // Read field on proc meshes
    PtrList<cloud> procClouds(procMeshes.size());
    PtrList<unallocatedIOPosition> procFields(procMeshes.size());
    forAll(procFields, proci)
    {
        const unallocatedFvMesh& procMesh = procMeshes[proci];

        Pout<< incrIndent;

        // Construct empty cloud
        procClouds.set
        (
            proci,
            new cloud
            (
                procMesh.thisDb(),
                "kinematicCloud"
            )
        );

        procFields.set
        (
            proci,
            new unallocatedIOPosition
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    procClouds[proci],
                    IOobject::MUST_READ,    //IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );

        Pout<< decrIndent;
    }

    unallocatedIOPosition particles
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    forAll(procFields, proci)
    {
        const unallocatedIOPosition& procCloud = procFields[proci];
        const labelList& cellMap = reconstructor.cellProcAddressing()[proci];

        forAllConstIter(typename IDLList<basicParticle>, procCloud, iter)
        {
            const basicParticle& p = iter();
            particles.append(new basicParticle(p, cellMap[p.cell()]));
        }
    }

    particles.writeData(os);

    Pout<< "** Reconstructed " << particles.size()
        << " into " << io.objectPath() << endl;
    return os.good();
}


bool Foam::passiveParticleStreamReconstructor::decompose
(
    const parUnallocatedFvFieldReconstructor& reconstructor,
    const unallocatedFvMesh& baseMesh,
    const IOobject& baseIO,

    const unallocatedFvMesh& thisMesh,
    const IOobject& thisIO,
    const bool,
    Ostream& os
) const
{
    Pout<< "*** LAGRANGIAN DEcomposing " << baseIO.objectPath() << endl;
    Pout<< "** LAGRANGIAN Decomposed " << baseIO.objectPath() << endl;
    return os.good();
}


bool Foam::passiveParticleStreamReconstructor::reconstruct
(
    const parUnallocatedFvFieldReconstructor& reconstructor,
    const regIOobject& thisIO,
    const bool,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    Pout<< "*** LAGRANGIAN reconstructin " << thisIO.objectPath() << endl;
    Pout<< "** LAGRANGIAN reconstruct " << thisIO.objectPath() << endl;
    return true;
}


// ************************************************************************* //
