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

#include "autoReconstructingFileOperation.H"
#include "uVolFields.H"
#include "unallocatedFvMesh.H"
#include "unallocatedVolMesh.H"
#include "unallocatedGenericFvPatchField.H"
#include "unallocatedGenericFvsPatchField.H"
#include "unallocatedFvFieldReconstructor.H"
#include "uFieldReconstructor.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::Ostream&
Foam::fileOperations::autoReconstructingFileOperation::
writeReconstructedFvVolumeField
(
    const fvMesh& mesh,
    const IOobject& io,
    Ostream& os
) const
{
    typedef GeometricField<Type, unallocatedFvPatchField, unallocatedVolMesh>
        GeoField;


    const uFieldReconstructor& reconstructor = uFieldReconstructor::New(mesh);

    const PtrList<unallocatedFvMesh>& procMeshes = reconstructor.procMeshes();

    // Read field on proc meshes
    PtrList<GeoField> procFields(procMeshes.size());
    forAll(procFields, proci)
    {
        const unallocatedFvMesh& procMesh = procMeshes[proci];
        procFields.set
        (
            proci,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    procMesh.time().timeName(),
                    procMesh.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMesh
            )
        );
    }

    // Fix filtering of empty nonuniform entries
    reconstructor.reconstructor().fixGenericNonuniform
    <
        GeoField,
        unallocatedGenericFvPatchField<Type>
    >(procFields);

    // Map local field onto baseMesh
    const unallocatedFvMesh& baseMesh = reconstructor.baseMesh();

    tmp<GeoField> tfld
    (
        reconstructor.reconstructor().reconstructFvVolumeField
        (
            IOobject
            (
                io.name(),
                baseMesh.time().timeName(),
                baseMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            procFields
        )
    );

    return os << tfld();
}


template<class Type>
Foam::Ostream&
Foam::fileOperations::autoReconstructingFileOperation::
writeReconstructedFvSurfaceField
(
    const fvMesh& mesh,
    const IOobject& io,
    Ostream& os
) const
{
    typedef GeometricField
    <
        Type,
        unallocatedFvsPatchField,
        unallocatedSurfaceMesh
    > GeoField;


    const uFieldReconstructor& reconstructor = uFieldReconstructor::New(mesh);

    const PtrList<unallocatedFvMesh>& procMeshes = reconstructor.procMeshes();

    // Read field on proc meshes
    PtrList<GeoField> procFields(procMeshes.size());
    forAll(procFields, proci)
    {
        const unallocatedFvMesh& procMesh = procMeshes[proci];
        procFields.set
        (
            proci,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    procMesh.time().timeName(),
                    procMesh.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMesh
            )
        );
    }

    // Fix filtering of empty nonuniform entries
    reconstructor.reconstructor().fixGenericNonuniform
    <
        GeoField,
        unallocatedGenericFvsPatchField<Type>
    >(procFields);

    // Map local field onto baseMesh
    const unallocatedFvMesh& baseMesh = reconstructor.baseMesh();

    tmp<GeoField> tfld
    (
        reconstructor.reconstructor().reconstructFvSurfaceField
        (
            IOobject
            (
                io.name(),
                baseMesh.time().timeName(),
                baseMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            procFields
        )
    );

    return os << tfld();
}


// ************************************************************************* //
