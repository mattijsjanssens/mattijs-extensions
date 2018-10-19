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

#include "unallocatedFvFieldReconstructor.H"
#include "fvFieldReconstructor.H"
#include "unallocatedGenericFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoField, class GenericType>
Foam::HashTable<Foam::wordPairList>
Foam::unallocatedFvFieldReconstructor::extractGenericTypes
(
    const GeoField& fld
) const
{
    // PatchField type to list of entry+type
    HashTable<wordPairList> patchFieldToTypes(2*fld.boundaryField().size());

    const typename GeoField::Boundary& bfld = fld.boundaryField();
    forAll(bfld, patchi)
    {
        if (isA<GenericType>(bfld[patchi]))
        {
            const GenericType& pfld =
                refCast<const GenericType>
                (
                    bfld[patchi]
                );

            genericPatchFieldBase::mergeEntries
            (
                pfld.actualTypeName(),
                pfld.entryTypes(),
                patchFieldToTypes
            );
        }
    }
    return patchFieldToTypes;
}


template<class GeoField, class GenericType>
void Foam::unallocatedFvFieldReconstructor::fixGenericNonuniform
(
    PtrList<GeoField>& procFields
) const
{
    // Extract field types of all genericFvPatchFields

    // PatchField type to list of entry+type
    HashTable<wordPairList> patchFieldToTypes
    (
        2*procFields[0].boundaryField().size()
    );

    forAll(procFields, proci)
    {
        HashTable<wordPairList> procTypes
        (
            extractGenericTypes<GeoField, GenericType>(procFields[proci])
        );

        forAllConstIter(HashTable<wordPairList>, procTypes, iter)
        {
            genericPatchFieldBase::mergeEntries
            (
                iter.key(),
                iter(),
                patchFieldToTypes
            );
        }
    }


    //DebugVar(patchFieldToTypes);


    // Add these to fields with unknown types
    forAll(procFields, proci)
    {
        typename GeoField::Boundary& bfld =
            procFields[proci].boundaryFieldRef();
        forAll(bfld, patchi)
        {
            if (isA<GenericType>(bfld[patchi]))
            {
                GenericType& pfld =
                    refCast<GenericType>
                    (
                        bfld[patchi]
                    );

                const wordPairList entryTypes(pfld.entryTypes());

                const wordPairList& allTypes =
                    patchFieldToTypes[pfld.actualTypeName()];

                forAll(allTypes, entryi)
                {
                    const Pair<word>& entry = allTypes[entryi];

                    if
                    (
                        genericPatchFieldBase::findEntry
                        (
                            entryTypes,
                            entry.first()
                        )
                     == -1
                    )
                    {
                        pfld.addEntry(entry.first(), entry.second());
                    }
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::unallocatedVolMesh>>
Foam::unallocatedFvFieldReconstructor::reconstructFvVolumeInternalField
(
    const IOobject& fieldIoObject,
    const PtrList<DimensionedField<Type, unallocatedVolMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nCells());

    forAll(procMeshes_, proci)
    {
        const DimensionedField<Type, unallocatedVolMesh>& procField =
            procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.field(),
            cellProcAddressing_[proci]
        );
    }

    return tmp<DimensionedField<Type, unallocatedVolMesh>>
    (
        new DimensionedField<Type, unallocatedVolMesh>
        (
            fieldIoObject,
            mesh_,
            procFields[0].dimensions(),
            internalField
        )
    );
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        Type,
        Foam::unallocatedFvPatchField,
        Foam::unallocatedVolMesh
    >
>
Foam::unallocatedFvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject,
    const PtrList
    <
        GeometricField
        <
            Type,
            unallocatedFvPatchField,
            unallocatedVolMesh
        >
    >& procFields
) const
{
    typedef GeometricField
    <
        Type,
        unallocatedFvPatchField,
        unallocatedVolMesh
    > GeoField;

    // Create the internalField
    Field<Type> internalField(mesh_.nCells());

    // Create the patch fields
    PtrList<unallocatedFvPatchField<Type>>
        patchFields(mesh_.boundary().size());

    forAll(procFields, proci)
    {
        const GeoField& procField = procFields[proci];

        if (procField.size() != procField.mesh().nCells())
        {
            FatalErrorInFunction<< "field size:" << procField.size()
                << " nCells:" << procField.mesh().nCells() << exit(FatalError);
        }
        if (procField.size() != cellProcAddressing_[proci].size())
        {
            FatalErrorInFunction<< "field size:" << procField.size()
                << " addressing:" << cellProcAddressing_[proci].size()
                << exit(FatalError);
        }

        //if (max(cellProcAddressing_[proci]) >= mesh_.nCells())
        //{
        //    FatalErrorInFunction<< "nTotalCells:" << mesh_.nCells()
        //        << " max addressing:" << max(cellProcAddressing_[proci])
        //        << exit(FatalError);
        //}


        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.primitiveField(),
            cellProcAddressing_[proci]
        );

         // Set the boundary patch values in the reconstructed field
         forAll(boundaryProcAddressing_[proci], patchi)
         {
             // Get patch index of the original patch
             const label curBPatch = boundaryProcAddressing_[proci][patchi];

            // Get addressing slice for this patch
            const labelList::subList cp =
                procField.mesh().boundary()[patchi].patchSlice
                (
                    faceProcAddressing_[proci]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                // Regular patch. Fast looping

                if (!patchFields(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        unallocatedFvPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, unallocatedVolMesh>::null(),
                            fvFieldReconstructor::fvPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                const label curPatchStart = mesh_.boundary()[curBPatch].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, facei)
                {
                    // Check
                    if (cp[facei] <= 0)
                    {
                        FatalErrorInFunction
                            << "Processor " << proci
                            << " patch "
                            << procField.mesh().boundary()[patchi].name()
                            << " face " << facei
                            << " originates from reversed face since "
                            << cp[facei]
                            << exit(FatalError);
                    }

                    // Subtract one to take into account offsets for
                    // face direction.
                    reverseAddressing[facei] = cp[facei] - 1 - curPatchStart;
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    reverseAddressing
                );
            }
            else
            {
                const Field<Type>& curProcPatch =
                    procField.boundaryField()[patchi];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, facei)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
                    label curF = cp[facei] - 1;

                    // Is the face on the boundary?
                    if (curF >= mesh_.nInternalFaces())
                    {
                        label curBPatch = mesh_.boundary().whichPatch(curF);

                        if (!patchFields(curBPatch))
                        {
                            patchFields.set
                            (
                                curBPatch,
                                unallocatedFvPatchField<Type>::New
                                (
                                    mesh_.boundary()[curBPatch].type(),
                                    mesh_.boundary()[curBPatch],
                                    DimensionedField
                                    <
                                        Type,
                                        unallocatedVolMesh
                                    >::null()
                                )
                            );
                        }

                        // add the face
                        label curPatchFace =
                            mesh_.boundary()[curBPatch].whichFace(curF);

                        patchFields[curBPatch][curPatchFace] =
                            curProcPatch[facei];
                    }
                }
            }
         }
    }


    //- Commented out - this is not needed when reconstructing 'normal' fields
    // forAll(mesh_.boundary(), patchi)
    // {
    //     // add empty patches
    //     if
    //     (
    //         isType<emptyFvPatch>(mesh_.boundary()[patchi])
    //      && !patchFields(patchi)
    //     )
    //     {
    //         patchFields.set
    //         (
    //             patchi,
    //             unallocatedFvPatchField<Type>::New
    //             (
    //                 emptyFvPatchField<Type>::typeName,
    //                 mesh_.boundary()[patchi],
    //                 DimensionedField<Type, unallocatedVolMesh>::null()
    //             )
    //         );
    //     }
    // }


    // Now construct from internalField and patchFields
    return tmp<GeoField>
    (
        new GeoField
        (
            fieldIoObject,
            mesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        Type,
        Foam::unallocatedFvsPatchField,
        Foam::unallocatedSurfaceMesh
    >
>
Foam::unallocatedFvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject,
    const PtrList
    <
        GeometricField
        <
            Type,
            unallocatedFvsPatchField,
            unallocatedSurfaceMesh
        >
    >& procFields
) const
{
    typedef GeometricField
    <
        Type,
        unallocatedFvsPatchField,
        unallocatedSurfaceMesh
    > GeoField;

    // Create the internalField
    Field<Type> internalField(mesh_.nInternalFaces());

    // Create the patch fields
    PtrList<unallocatedFvsPatchField<Type>> patchFields
    (
        mesh_.boundary().size()
    );

    forAll(procFields, proci)
    {
        const GeoField& procField = procFields[proci];

        if (procField.size() != procField.mesh().nInternalFaces())
        {
            FatalErrorInFunction<< "field size:" << procField.size()
                << " nInternalFaces:" << procField.mesh().nInternalFaces()
                << exit(FatalError);
        }

        // Set the face values in the reconstructed field

        // It is necessary to create a copy of the addressing array to
        // take care of the face direction offset trick.
        //
        {
            const labelList& faceMap = faceProcAddressing_[proci];

            // Correctly oriented copy of internal field
            Field<Type> procInternalField(procField.primitiveField());

            // Addressing into original field
            labelList curAddr(procInternalField.size());

            forAll(procInternalField, addrI)
            {
                curAddr[addrI] = mag(faceMap[addrI])-1;
                if (faceMap[addrI] < 0)
                {
                    procInternalField[addrI] = -procInternalField[addrI];
                }
            }

            // Map
            internalField.rmap(procInternalField, curAddr);
        }

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[proci], patchi)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[proci][patchi];

            // Get addressing slice for this patch
            const labelList::subList cp =
                procField.mesh().boundary()[patchi].patchSlice
                (
                    faceProcAddressing_[proci]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                // Regular patch. Fast looping

                if (!patchFields(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        unallocatedFvsPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            DimensionedField
                            <
                                Type,
                                unallocatedSurfaceMesh
                            >::null(),
                            fvFieldReconstructor::fvPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                const label curPatchStart = mesh_.boundary()[curBPatch].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, facei)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
                    reverseAddressing[facei] = cp[facei] - 1 - curPatchStart;
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    reverseAddressing
                );
            }
            else
            {
                const Field<Type>& curProcPatch =
                    procField.boundaryField()[patchi];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, facei)
                {
                    label curF = cp[facei] - 1;

                    // Is the face turned the right side round
                    if (curF >= 0)
                    {
                        // Is the face on the boundary?
                        if (curF >= mesh_.nInternalFaces())
                        {
                            label curBPatch = mesh_.boundary().whichPatch(curF);

                            if (!patchFields(curBPatch))
                            {
                                patchFields.set
                                (
                                    curBPatch,
                                    unallocatedFvsPatchField<Type>::New
                                    (
                                        mesh_.boundary()[curBPatch].type(),
                                        mesh_.boundary()[curBPatch],
                                        DimensionedField
                                        <
                                            Type,
                                            unallocatedSurfaceMesh
                                        >::null()
                                    )
                                );
                            }

                            // add the face
                            label curPatchFace =
                                mesh_.boundary()[curBPatch].whichFace(curF);

                            patchFields[curBPatch][curPatchFace] =
                                curProcPatch[facei];
                        }
                        else
                        {
                            // Internal face
                            internalField[curF] = curProcPatch[facei];
                        }
                    }
                }
            }
        }
    }

//     forAll(mesh_.boundary(), patchi)
//     {
//         // add empty patches
//         if
//         (
//             isType<emptyFvsPatch>(mesh_.boundary()[patchi])
//          && !patchFields(patchi)
//         )
//         {
//             patchFields.set
//             (
//                 patchi,
//                 unallocatedFvsPatchField<Type>::New
//                 (
//                     emptyunallocatedFvsPatchField<Type>::typeName,
//                     mesh_.boundary()[patchi],
//                     DimensionedField<Type, unallocatedSurfaceMesh>::null()
//                 )
//             );
//         }
//     }


    // Now construct and write the field
    // setting the internalField and patchFields
    return tmp<GeoField>
    (
        new GeoField
        (
            fieldIoObject,
            mesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


// ************************************************************************* //
