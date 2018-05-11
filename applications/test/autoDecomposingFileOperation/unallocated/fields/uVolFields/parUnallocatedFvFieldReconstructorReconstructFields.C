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
#include "Time.H"
#include "PtrList.H"
#include "fvPatchFields.H"
#include "emptyFvPatchField.H"
#include "processorFvPatchField.H"
#include "IOobjectList.H"
#include "mapDistributePolyMesh.H"
#include "processorFvPatch.H"

#include "directFvPatchFieldMapper.H"
#include "distributedUnallocatedDirectFieldMapper.H"
#include "distributedUnallocatedDirectFvPatchFieldMapper.H"

#include "fvsPatchFields.H"
#include "unallocatedFvMesh.H"
#include "unallocatedFvBoundaryMesh.H"
#include "unallocatedSurfaceMesh.H"

#include "distributedDirectFieldMapper.H"
#include "distributedDirectFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type>
// Foam::tmp<Foam::DimensionedField<Type, Foam::unallocatedVolMesh> >
// Foam::parUnallocatedFvFieldReconstructor::reconstructFvVolumeInternalField
// (
//     const DimensionedField<Type, volMesh>& fld
// ) const
// {
//     distributedUnallocatedDirectFieldMapper mapper
//     (
//         labelUList::null(),
//         distMap_.cellMap()
//     );
//
//     Field<Type> internalField(fld, mapper);
//
//     // Construct a volField
//     IOobject baseIO
//     (
//         fld.name(),
//         baseMesh_.time().timeName(),
//         fld.local(),
//         baseMesh_,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     );
//
//     return tmp<DimensionedField<Type, unallocatedVolMesh> >
//     (
//         new DimensionedField<Type, unallocatedVolMesh>
//         (
//             baseIO,
//             baseMesh_,
//             fld.dimensions(),
//             internalField
//         )
//     );
// }


template<class GeoField>
Foam::tmp<GeoField>
Foam::parUnallocatedFvFieldReconstructor::reconstructFvVolumeField
(
    const GeoField& fld
) const
{
    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap()
    );

    Field<typename GeoField::value_type> internalField
    (
        fld.internalField(),
        mapper
    );

    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on mesh, not baseMesh

    PtrList<typename GeoField::Patch> patchFields
    (
        fld.mesh().boundary().size()
    );

    const typename GeoField::Boundary& bfld = fld.boundaryField();

    forAll(bfld, patchI)
    {
        DebugVar(bfld[patchI].patch().name());

        if (patchFaceMaps_.set(patchI))
        {
            Pout<< "** mapping patch " << bfld[patchI].patch().name()
                << " patchField type:" << bfld[patchI].type()
                << " of field " << fld.name() << endl;

            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchI]
            );

            // Map into local copy
            patchFields[patchI].autoMap(mapper);
        }
    }


    PtrList<typename GeoField::Patch> basePatchFields
    (
        baseMesh_.boundary().size()
    );

    // Clone the patchFields onto the base patches. This is just to reset
    // the reference to the patch, size and content stay the same.
    forAll(patchFields, patchI)
    {
        if (patchFields.set(patchI))
        {
            const fvPatch& basePatch = baseMesh_.boundary()[patchI];

            const typename GeoField::Patch& pfld = patchFields[patchI];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            basePatchFields.set
            (
                patchI,
                GeoField::Patch::New
                (
                    pfld,
                    basePatch,
                    GeoField::Internal::null(),
                    dummyMapper
                )
            );
        }
    }

    // // Add some empty patches on remaining patches (tbd.probably processor
    // // patches)
    // forAll(basePatchFields, patchI)
    // {
    //     if (patchI >= patchFields.size() || !patchFields.set(patchI))
    //     {
    //         basePatchFields.set
    //         (
    //             patchI,
    //             GeoField::Patch::New
    //             (
    //                 baseMesh_.boundary()[patchI].type(),
    //                 baseMesh_.boundary()[patchI],
    //                 GeoField::Internal::null()
    //             )
    //         );
    //     }
    // }

    // Construct a GeoField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_.thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    return tmp<GeoField>
    (
        new GeoField
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::parUnallocatedFvFieldReconstructor::reconstructFvSurfaceField
(
    const GeoField& fld
) const
{
    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Note: faceMap() is for all faces, not just internal ones so
    //       just extend field before mapping. Easier than rewriting
    //       faceMap to be internal-to-internal

    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.faceMap()
    );

    Field<typename GeoField::value_type> intFld(fld.internalField());
    // Pad to full length
    intFld.setSize
    (
        distMap_.faceMap().constructSize(),
        pTraits<typename GeoField::value_type>::zero
    );
    // Do mapping
    Field<typename GeoField::value_type> internalField(intFld, mapper);
    // Shrink back to wanted size
    internalField.setSize(unallocatedSurfaceMesh::size(baseMesh_));


    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on mesh, not baseMesh

    PtrList<typename GeoField::Patch> patchFields
    (
        fld.mesh().boundary().size()
    );

    const typename GeoField::Boundary& bfld = fld.boundaryField();

    forAll(bfld, patchI)
    {
        DebugVar(bfld[patchI].patch().name());

        if (patchFaceMaps_.set(patchI))
        {
            Pout<< "** mapping patch " << bfld[patchI].patch().name()
                << " patchField type:" << bfld[patchI].type()
                << " of field " << fld.name() << endl;

            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchI]
            );

            // Map into local copy
            patchFields[patchI].autoMap(mapper);
        }
    }


    PtrList<typename GeoField::Patch> basePatchFields
    (
        baseMesh_.boundary().size()
    );

    // Clone the patchFields onto the base patches. This is just to reset
    // the reference to the patch, size and content stay the same.
    forAll(patchFields, patchI)
    {
        if (patchFields.set(patchI))
        {
            const fvPatch& basePatch = baseMesh_.boundary()[patchI];

            const typename GeoField::Patch& pfld = patchFields[patchI];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            basePatchFields.set
            (
                patchI,
                GeoField::Patch::New
                (
                    pfld,
                    basePatch,
                    GeoField::Internal::null(),
                    dummyMapper
                )
            );
        }
    }

    // // Add some empty patches on remaining patches (tbd.probably processor
    // // patches)
    // forAll(basePatchFields, patchI)
    // {
    //     if (patchI >= patchFields.size() || !patchFields.set(patchI))
    //     {
    //         basePatchFields.set
    //         (
    //             patchI,
    //             GeoField::Patch::New
    //             (
    //                 baseMesh_.boundary()[patchI].type(),
    //                 baseMesh_.boundary()[patchI],
    //                 GeoField::Internal::null()
    //             )
    //         );
    //     }
    // }

    // Construct a GeoField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_.thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    return tmp<GeoField>
    (
        new GeoField
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );
}


//XXXXXXXXXX
template<class GeoField>
Foam::tmp<GeoField>
Foam::parUnallocatedFvFieldReconstructor::decomposeFvVolumeField
(
    const GeoField& fld
) const
{
    if (fld.size() != baseMesh_.nCells())
    {
        FatalErrorInFunction<< "Size:" << fld.size()
            << " base size:" << baseMesh_.nCells() << exit(FatalError);
    }


    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Create reverse mapper
    distributedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap(),
        procMesh_.nCells()              // Construct size
    );

    Field<typename GeoField::value_type> internalField
    (
        fld.internalField(),
        mapper
    );

    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on baseMesh, not procMesh

    PtrList<typename GeoField::Patch> patchFields
    (
        fld.mesh().boundary().size()
    );

    const typename GeoField::Boundary& bfld = fld.boundaryField();

    forAll(bfld, patchI)
    {
        DebugVar(bfld[patchI].patch().name());

        if (patchFaceMaps_.set(patchI))
        {
            Pout<< "** mapping patch " << bfld[patchI].patch().name()
                << " patchField type:" << bfld[patchI].type()
                << " of field " << fld.name() << endl;

            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            Pout<< "** patchField size:" << patchFields[patchI].size()
                << endl;

            distributedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchI],
                procMesh_.boundary()[patchI].size()
            );

            // Map into local copy
            patchFields[patchI].autoMap(mapper);
        }
    }


    PtrList<typename GeoField::Patch> procPatchFields
    (
        procMesh_.boundary().size()
    );

    // Clone the patchFields onto the proc patches. This is just to reset
    // the reference to the patch, size and content stay the same.
    forAll(patchFields, patchI)
    {
        if (patchFields.set(patchI))
        {
            const fvPatch& procPatch = procMesh_.boundary()[patchI];

            const typename GeoField::Patch& pfld = patchFields[patchI];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            procPatchFields.set
            (
                patchI,
                GeoField::Patch::New
                (
                    pfld,
                    procPatch,
                    GeoField::Internal::null(),
                    dummyMapper
                )
            );
        }
    }

    // Add some empty patches on remaining patches (tbd.probably processor
    // patches)
    forAll(procPatchFields, patchI)
    {
        if (!procPatchFields.set(patchI))
        {
            const fvPatch& procPatch = procMesh_.boundary()[patchI];

Pout<< "** synthesising field of type " << procPatch.type()
    << " on " << procPatch.name() << endl;

            word patchType(procPatch.type());
//             if (isA<unallocatedGenericFvPatch>(procPatch))
//             {
//                 patchType = refCast
//                 <
//                     const unallocatedGenericFvPatch
//                 >(procPatch).actualTypeName();
//             }
DebugVar(patchType);


            procPatchFields.set
            (
                patchI,
                GeoField::Patch::New
                (
                    patchType,
                    procPatch,
                    GeoField::Internal::null()
                )
            );
        }
    }

    // Construct a GeoField
    IOobject procIO
    (
        fld.name(),
        procMesh_.time().timeName(),
        fld.local(),
        procMesh_.thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    return tmp<GeoField>
    (
        new GeoField
        (
            procIO,
            procMesh_,
            fld.dimensions(),
            internalField,
            procPatchFields
        )
    );
}
//XXXXXXXXXX

// ************************************************************************* //
