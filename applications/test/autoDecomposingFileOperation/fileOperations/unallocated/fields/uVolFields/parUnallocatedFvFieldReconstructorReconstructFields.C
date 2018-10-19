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
#include "SubField.H"

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

    forAll(baseMesh_.boundary(), patchI)
    {
        if (patchReconFaceMaps_.set(patchI))
        {
            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchReconFaceMaps_[patchI]
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
    forAll(baseMesh_.boundary(), patchI)
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

    // Construct a GeoField
    return tmp<GeoField>
    (
        new GeoField
        (
            IOobject
            (
                fld.name(),
                baseMesh_.time().timeName(),
                fld.local(),
                baseMesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::parUnallocatedFvFieldReconstructor::flatten
(
    Field<Type>& flatFld,
    const GeometricField<Type, PatchField, GeoMesh>& fld
)
{
    const typename GeometricField<Type, PatchField, GeoMesh>::Boundary&
        bfld = fld.boundaryField();

    SubField<Type>
    (
        flatFld,
        fld.internalField().size()
    ) = fld.internalField();


    forAll(bfld, patchI)
    {
        SubField<Type>
        (
            flatFld,
            bfld[patchI].size(),
            bfld[patchI].patch().start()
        ) = bfld[patchI];
    }
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::parUnallocatedFvFieldReconstructor::reconstructFvSurfaceField
(
    const GeoField& fld
) const
{
    // surface fields need to map the processor fields back to internal
    // faces. Two choices:
    // - flatten all boundary values and use full map
    // - create map to override internal from processor faces
    // For now use flattening

    const typename GeoField::Boundary& bfld = fld.boundaryField();


    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.faceMap()
    );

    // Create internal values + boundary values
    Field<typename GeoField::value_type> flatFld
    (
        distMap_.faceMap().constructSize()
    );
    flatten(flatFld, fld);

    // Do actual mapping
    Field<typename GeoField::value_type> internalField(flatFld, mapper);
    // Shrink back to wanted size
    internalField.setSize(baseMesh_.nInternalFaces());


    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on mesh, not baseMesh

    PtrList<typename GeoField::Patch> patchFields
    (
        fld.mesh().boundary().size()
    );

    forAll(baseMesh_.boundary(), patchI)
    {
        if (patchReconFaceMaps_.set(patchI))
        {
            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            const distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchReconFaceMaps_[patchI]
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
    forAll(baseMesh_.boundary(), patchI)
    {
        if (patchFields.set(patchI))
        {
            const fvPatch& basePatch = baseMesh_.boundary()[patchI];
            const typename GeoField::Patch& pfld = patchFields[patchI];

            const labelList dummyMap(identity(pfld.size()));
            const directFvPatchFieldMapper dummyMapper(dummyMap);

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

    // Construct a GeoField
    return tmp<GeoField>
    (
        new GeoField
        (
            IOobject
            (
                fld.name(),
                baseMesh_.time().timeName(),
                fld.local(),
                baseMesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );
}


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

    // This is a bit overkill - all processors have the full undecomposed
    // field so there is no need for remote access ... However at some
    // point maybe only master reads.

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
        if (patchDecompFaceMaps_.set(patchI))
        {
            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            const distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchDecompFaceMaps_[patchI]
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

            const labelList dummyMap(identity(pfld.size()));
            const directFvPatchFieldMapper dummyMapper(dummyMap);

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


    // Processor patches
    // -----------------
    // Add processor patches on remaining patches. The processor patches
    // are special in that their value is mapped from (remote) internal values.
    // So:
    // - swap internal values next to processor patches
    // - construct processor patch fields with these remote values

    typedef typename GeoField::value_type Type;

    List<Field<Type>> localValues(Pstream::nProcs());
    List<Field<Type>> remoteValues(Pstream::nProcs());

    forAll(procPatchFields, patchI)
    {
        if (!procPatchFields.set(patchI))
        {
            const unallocatedGenericFvPatch& pp = procMesh_.boundary()[patchI];
            if (pp.type() == processorFvPatch::typeName)
            {
                // Use the dictionary to lookup nbr
                label nbrProci = readLabel(pp.dict().lookup("neighbProcNo"));
                localValues[nbrProci] = Field<Type>
                (
                    internalField,
                    pp.faceCells()
                );
            }
        }
    }

    labelList localSizes(Pstream::nProcs(), 0);
    forAll(localValues, proci)
    {
        localSizes[proci] = localValues[proci].size();
    }

    Pstream::exchange<Field<Type>, Type>
    (
        localValues,
        localSizes,
        remoteValues
    );

    forAll(procPatchFields, patchI)
    {
        if (!procPatchFields.set(patchI))
        {
            const unallocatedGenericFvPatch& pp = procMesh_.boundary()[patchI];
            if (pp.type() == processorFvPatch::typeName)
            {
                label nbrProci = readLabel(pp.dict().lookup("neighbProcNo"));

                procPatchFields.set
                (
                    patchI,
                    GeoField::Patch::New
                    (
                        pp.type(),
                        pp,
                        GeoField::Internal::null()
                    )
                );
                procPatchFields[patchI] = remoteValues[nbrProci];
            }
        }
    }

    // Construct a GeoField
    return tmp<GeoField>
    (
        new GeoField
        (
            IOobject
            (
                fld.name(),
                procMesh_.time().timeName(),
                fld.local(),
                procMesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            fld.dimensions(),
            internalField,
            procPatchFields
        )
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::parUnallocatedFvFieldReconstructor::decomposeFvSurfaceField
(
    const GeometricField<Type, PatchField, GeoMesh>& fld
) const
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    if (fld.size() != GeoMesh::size(baseMesh_))
    {
        FatalErrorInFunction<< "Size:" << fld.size()
            << " base size:" << GeoMesh::size(baseMesh_) << exit(FatalError);
    }

    const typename GeoField::Boundary& bfld = fld.boundaryField();


    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // See comment above about using map over all faces (faceMap) v.s.
    // creating a 'proper' map over the internal faces only

    // Create reverse mapper
    distributedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.faceMap(),
        procMesh_.nFaces()              // Construct size
    );

    // Create internal values + boundary values
    Field<typename GeoField::value_type> flatFld
    (
        distMap_.faceMap().constructSize()
    );
    flatten(flatFld, fld);

    // Do actual mapping
    // Now we have a value for all faces (internal or boundary) on the proc mesh
    Field<typename GeoField::value_type> internalField(flatFld, mapper);


    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on baseMesh, not procMesh

    PtrList<typename GeoField::Patch> patchFields
    (
        fld.mesh().boundary().size()
    );

    forAll(bfld, patchI)
    {
        if (patchDecompFaceMaps_.set(patchI))
        {
            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            const distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchDecompFaceMaps_[patchI]
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

            const labelList dummyMap(identity(pfld.size()));
            const directFvPatchFieldMapper dummyMapper(dummyMap);

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

    // Add processorPatchFields on remaining patches
    forAll(procPatchFields, patchI)
    {
        if (!procPatchFields.set(patchI))
        {
            const unallocatedGenericFvPatch& pp = procMesh_.boundary()[patchI];
            if (pp.type() == processorFvPatch::typeName)
            {
                procPatchFields.set
                (
                    patchI,
                    GeoField::Patch::New
                    (
                        pp.type(),
                        pp,
                        GeoField::Internal::null()
                    )
                );

                procPatchFields[patchI] =
                    SubField<typename GeoField::value_type>
                    (
                        internalField,
                        pp.size(),
                        pp.start()
                    );
            }
        }
    }

    // Shrink field back to internal faces only now we've used the processor
    // values
    internalField.setSize(GeoMesh::size(procMesh_));

    // Construct a GeoField
    return tmp<GeoField>
    (
        new GeoField
        (
            IOobject
            (
                fld.name(),
                procMesh_.time().timeName(),
                fld.local(),
                procMesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            fld.dimensions(),
            internalField,
            procPatchFields
        )
    );
}


// ************************************************************************* //
