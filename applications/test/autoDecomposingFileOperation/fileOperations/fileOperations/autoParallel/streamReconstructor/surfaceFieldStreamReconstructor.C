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

#include "surfaceFieldStreamReconstructor.H"
#include "uSurfaceFields.H"
#include "unallocatedFvMesh.H"
#include "unallocatedSurfaceMesh.H"
#include "unallocatedGenericFvsPatchField.H"
#include "unallocatedFvFieldReconstructor.H"
#include "uFieldReconstructor.H"
#include "parUnallocatedFvFieldReconstructor.H"
#include "unallocatedFvMeshObject.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::surfaceFieldStreamReconstructor<Type>::reconstruct
(
    const IOobject& io,
    const bool,
    Ostream& os
) const
{
    typedef GeometricField
    <
        Type,
        unallocatedFvsPatchField,
        unallocatedSurfaceMesh
    > GeoField;

    const uFieldReconstructor& reconstructor =
        uFieldReconstructor::New(io.db());

    const PtrList<unallocatedFvMesh>& procMeshes = reconstructor.procMeshes();

    Info<< "Reconstructing " << io.objectPath() << endl;

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
                    io.instance(),
                    io.local(),
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
                io.instance(),
                io.local(),
                baseMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            procFields
        )
    );

    Pout<< incrIndent;
    os << tfld();
    Pout<< decrIndent;

    return os.good();
}


template<class Type>
bool Foam::surfaceFieldStreamReconstructor<Type>::decompose
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
    typedef GeometricField
    <
        Type,
        unallocatedFvsPatchField,
        unallocatedSurfaceMesh
    > GeoField;

    // Read base field
    Info<< "Reading " << baseIO.objectPath() << endl;
    const GeoField baseFld(baseIO, baseMesh);

    // Decompose
    tmp<GeoField> tfld
    (
        reconstructor.decomposeFvSurfaceField(baseFld)
    );

    // Stream
    Pout<< incrIndent;
    os << tfld();
    Pout<< decrIndent;

    return os.good();
}


template<class Type>
bool Foam::surfaceFieldStreamReconstructor<Type>::reconstruct
(
    const parUnallocatedFvFieldReconstructor& reconstructor,
    const regIOobject& thisIO,
    const bool,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> GeoField;

    const GeoField& fld = dynamic_cast<const GeoField&>(thisIO);
    const unallocatedFvMesh& uMesh = unallocatedFvMeshObject::New(fld.mesh());

    // We have a slight problem: fld is real field on real mesh (fvMesh),
    // but the reconstructor is on an unallocated mesh. We cannot currently
    // reconstruct from a real mesh onto an unallocated mesh.
    // As a workaround convert the fld to its unallocated version

    OStringStream os(IOstream::BINARY);
    os  << fld;
    IStringStream is(os.str(), IOstream::BINARY);

    typedef GeometricField
    <
        Type,
        unallocatedFvsPatchField,
        unallocatedSurfaceMesh
    > UGeoField;
    const UGeoField uProcFld(fld, uMesh, dictionary(is));

    // Map local field onto baseMesh
    tmp<UGeoField> tfld
    (
        reconstructor.reconstructFvSurfaceField(uProcFld)
    );

    // Write master field to parent
    bool state = true;
    if (Pstream::master())
    {
        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;
        {
            Info<< "Writing " << tfld().objectPath() << endl;
            Pout<< incrIndent;
            state = tfld().writeObject(fmt, ver, cmp, true);
            Pout<< decrIndent;
        }
        Pstream::parRun() = oldParRun;
    }

    return state;
}


// ************************************************************************* //
