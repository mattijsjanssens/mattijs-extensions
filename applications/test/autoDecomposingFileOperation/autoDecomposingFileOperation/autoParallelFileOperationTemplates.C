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

#include "parUnallocatedFvFieldReconstructor.H"
#include "uVolFields.H"
#include "unallocatedFvMeshObject.H"
#include "uSurfaceFields.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoField>
bool Foam::fileOperations::autoParallelFileOperation::reconstructAndWrite
(
    const unallocatedFvMesh& uMesh,
    const GeoField& procFld,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    if (!reconstructorPtr_.valid())
    {
        // Read base mesh
        const unallocatedFvMesh& baseM = baseMesh(procFld.time());

        Info<< "Creating reconstructor" << nl << endl;
        reconstructorPtr_ = new parUnallocatedFvFieldReconstructor
        (
            baseM,
            uMesh,
            distMapPtr_()
        );
    }
    const parUnallocatedFvFieldReconstructor& reconstructor =
        reconstructorPtr_();

    // Map local field onto baseMesh
    tmp<GeoField> tfld
    (
        reconstructor.reconstructFvVolumeField(procFld)
    );

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


template<class GeoField>
bool Foam::fileOperations::autoParallelFileOperation::reconstructAndWrite
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    Pout<< "** reconstructing:" << io.objectPath() << endl;
    const GeoField& fld = dynamic_cast<const GeoField&>(io);

    const unallocatedFvMesh& uMesh = unallocatedFvMeshObject::New(fld.mesh());

    // We have a slight problem: fld is real field on real mesh (fvMesh),
    // but the reconstructor is on an unallocated mesh. We cannot currently
    // reconstruct from a real mesh onto an unallocated mesh.
    // As a workaround convert the fld to its unallocated version

    OStringStream os(IOstream::BINARY);
    os  << fld;
    IStringStream is(os.str(), IOstream::BINARY);

    GeometricField
    <
        typename GeoField::value_type,
        unallocatedFvPatchField,
        unallocatedVolMesh
    >
    uProcFld(fld, uMesh, dictionary(is));

    return reconstructAndWrite(uMesh, uProcFld, fmt, ver, cmp);
}


template<class GeoField>
bool Foam::fileOperations::autoParallelFileOperation::reconstructAndWrite2
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    Pout<< "** reconstructing:" << io.objectPath() << endl;
    const GeoField& procFld = dynamic_cast<const GeoField&>(io);

    const unallocatedFvMesh& uMesh =
        unallocatedFvMeshObject::New(procFld.mesh());

    // We have a slight problem: fld is real field on real mesh (fvMesh),
    // but the reconstructor is on an unallocated mesh. We cannot currently
    // reconstruct from a real mesh onto an unallocated mesh.
    // As a workaround convert the fld to its unallocated version

    typedef GeometricField
    <
        typename GeoField::value_type,
        unallocatedFvsPatchField,
        unallocatedSurfaceMesh
    > uGeoField;


    OStringStream os(IOstream::BINARY);
    os  << procFld;
    IStringStream is(os.str(), IOstream::BINARY);
    uGeoField uProcFld(procFld, uMesh, dictionary(is));


    if (!reconstructorPtr_.valid())
    {
        // Read base mesh
        const unallocatedFvMesh& baseM = baseMesh(uProcFld.time());

        Info<< "Creating reconstructor" << nl << endl;
        reconstructorPtr_ = new parUnallocatedFvFieldReconstructor
        (
            baseM,
            uMesh,
            distMapPtr_()
        );
    }
    const parUnallocatedFvFieldReconstructor& reconstructor =
        reconstructorPtr_();

    // Map local field onto baseMesh
    tmp<uGeoField> tfld
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
            state = tfld().writeObject(fmt, ver, cmp, true);
        }
        Pstream::parRun() = oldParRun;
    }

    return state;
}


// template<class GeoField>
// bool Foam::fileOperations::autoParallelFileOperation::decomposeAndWrite
// (
//     const IOobject& procIO,
//     const IOobject& parentIO,
//     const word& type,
//     Ostream& os
// ) const
// {
//     if (type == GeoField::typeName)
//     {
//         Info<< "Loading field " << parentIO.name() << endl;
//
//         GeometricField
//         <
//             typename GeoField::value_type,
//             unallocatedFvPatchField,
//             unallocatedVolMesh
//         > baseFld
//         (
//             parentIO,
// XXXXXX
//
//             IOobject
//             (
//                 "p",
//                 mesh.time().timeName(),
//                 mesh.thisDb(),
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE,
//                 false
//             ),
//             mesh
//         );
//
//
//         const fvMesh& undecomposedMesh = baseMesh(parentIO.time());
//         GeoField parentFld(parentIO, undecomposedMesh);
//
//         const fvFieldDecomposer& volDecomposer(decomposer(procIO));
//
//         // Decompose field and transfer to stream
//         tmp<GeoField> fld(volDecomposer.decomposeField(parentFld));
//         fld().writeData(os);
//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }

// ************************************************************************* //
