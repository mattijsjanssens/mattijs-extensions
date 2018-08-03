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

#include "volFieldStreamReconstructor.H"
#include "uVolFields.H"
#include "unallocatedFvMesh.H"
#include "unallocatedVolMesh.H"
#include "unallocatedGenericFvPatchField.H"
#include "unallocatedFvFieldReconstructor.H"
#include "uFieldReconstructor.H"
#include "parUnallocatedFvFieldReconstructor.H"
#include "unallocatedFvMeshObject.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::volFieldStreamReconstructor<Type>::reconstruct
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

    os << tfld();

    return os.good();
}


template<class Type>
bool Foam::volFieldStreamReconstructor<Type>::decompose
(
    const parUnallocatedFvFieldReconstructor& reconstructor,
    const unallocatedFvMesh& baseMesh,
    const IOobject& baseIO,

    const unallocatedFvMesh& thisMesh,
    const IOobject& thisIO,
    Ostream& os
) const
{
    typedef GeometricField<Type, unallocatedFvPatchField, unallocatedVolMesh>
        GeoField;

    // Read base field
    const GeoField baseFld(baseIO, baseMesh);

    // Decompose
    tmp<GeoField> tfld
    (
        reconstructor.decomposeFvVolumeField(baseFld)
    );

    // Stream
    os << tfld();

    return os.good();
}


template<class Type>
bool Foam::volFieldStreamReconstructor<Type>::reconstruct
(
    const parUnallocatedFvFieldReconstructor& reconstructor,
    const regIOobject& thisIO,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    const GeoField& fld = dynamic_cast<const GeoField&>(thisIO);
    const unallocatedFvMesh& uMesh = unallocatedFvMeshObject::New(fld.mesh());

    // We have a slight problem: fld is real field on real mesh (fvMesh),
    // but the reconstructor is on an unallocated mesh. We cannot currently
    // reconstruct from a real mesh onto an unallocated mesh.
    // As a workaround convert the fld to its unallocated version

    OStringStream os(IOstream::BINARY);
    os  << fld;
    IStringStream is(os.str(), IOstream::BINARY);

    typedef GeometricField<Type, unallocatedFvPatchField, unallocatedVolMesh>
        UGeoField;
    const UGeoField uProcFld(fld, uMesh, dictionary(is));

    // Map local field onto baseMesh
    tmp<UGeoField> tfld(reconstructor.reconstructFvVolumeField(uProcFld));

    // Write master field to parent
    bool state = true;
    if (Pstream::master())
    {
        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;
        {
            //Pout<< "**Writign " << tfld().objectPath() << endl;
            //state = tfld().write();
            state = tfld().writeObject(fmt, ver, cmp, true);
        }
        Pstream::parRun() = oldParRun;
    }

    return state;
}


// ************************************************************************* //
