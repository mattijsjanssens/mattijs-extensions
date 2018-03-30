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
#include "Time.H"
#include "PtrList.H"
#include "fvPatchFields.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchField.H"
#include "emptyFvsPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::unallocatedFvMesh>>
Foam::unallocatedFvFieldReconstructor::reconstructFvVolumeInternalField
(
    const IOobject& fieldIoObject,
    const PtrList<DimensionedField<Type, unallocatedFvMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nCells());

    forAll(procMeshes_, proci)
    {
        const DimensionedField<Type, unallocatedFvMesh>& procField =
            procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.field(),
            cellProcAddressing_[proci]
        );
    }

    return tmp<DimensionedField<Type, unallocatedFvMesh>>
    (
        new DimensionedField<Type, unallocatedFvMesh>
        (
            fieldIoObject,
            mesh_,
            procFields[0].dimensions(),
            internalField
        )
    );
}


//template<class Type>
//Foam::tmp
//<Foam::GeometricField<Type,
//Foam::unallocatedFvPatchField, Foam::unallocatedFvMesh>>
//Foam::unallocatedFvFieldReconstructor::reconstructFvVolumeField
//(
//    const IOobject& fieldIoObject,
//    const PtrList
//    <
//        GeometricField
//        <
//            Type,
//            unallocatedFvPatchField,
//            unallocatedFvMesh
//        >
//    >& procFields
//) const
//{
//    // Create the internalField
//    Field<Type> internalField(mesh_.nCells());
//
//    // Create the patch fields
//    PtrList<unallocatedFvPatchField<Type>>
//        patchFields(mesh_.boundary().size());
//
//    forAll(procFields, proci)
//    {
//        const GeometricField
//         <Type, unallocatedFvPatchField, unallocatedFvMesh>& procField =
//            procFields[proci];
//
//        // Set the cell values in the reconstructed field
//        internalField.rmap
//        (
//            procField.primitiveField(),
//            cellProcAddressing_[proci]
//        );
//
//        // Set the boundary patch values in the reconstructed field
//        forAll(boundaryProcAddressing_[proci], patchi)
//        {
//            // Get patch index of the original patch
//            const label curBPatch = boundaryProcAddressing_[proci][patchi];
//
//            // Get addressing slice for this patch
//            const labelList::subList cp =
//                procField.mesh().boundary()[patchi].patchSlice
//                (
//                    faceProcAddressing_[proci]
//                );
//
//            // check if the boundary patch is not a processor patch
//            if (curBPatch >= 0)
//            {
//                // Regular patch. Fast looping
//
//                if (!patchFields(curBPatch))
//                {
//                    patchFields.set
//                    (
//                        curBPatch,
//                        unallocatedFvPatchField<Type>::New
//                        (
//                            procField.boundaryField()[patchi],
//                            mesh_.boundary()[curBPatch],
//                            DimensionedField<Type, unallocatedFvMesh>::null(),
//                            unallocatedFvPatchFieldReconstructor
//                            (
//                                mesh_.boundary()[curBPatch].size()
//                            )
//                        )
//                    );
//                }
//
//                const label curPatchStart =
//                    mesh_.boundaryMesh()[curBPatch].start();
//
//                labelList reverseAddressing(cp.size());
//
//                forAll(cp, facei)
//                {
//                    // Check
//                    if (cp[facei] <= 0)
//                    {
//                        FatalErrorInFunction
//                            << "Processor " << proci
//                            << " patch "
//                            << procField.mesh().boundary()[patchi].name()
//                            << " face " << facei
//                            << " originates from reversed face since "
//                            << cp[facei]
//                            << exit(FatalError);
//                    }
//
//                    // Subtract one to take into account offsets for
//                    // face direction.
//                    reverseAddressing[facei] = cp[facei] - 1 - curPatchStart;
//                }
//
//
//                patchFields[curBPatch].rmap
//                (
//                    procField.boundaryField()[patchi],
//                    reverseAddressing
//                );
//            }
//            else
//            {
//                const Field<Type>& curProcPatch =
//                    procField.boundaryField()[patchi];
//
//                // In processor patches, there's a mix of internal faces (some
//                // of them turned) and possible cyclics. Slow loop
//                forAll(cp, facei)
//                {
//                    // Subtract one to take into account offsets for
//                    // face direction.
//                    label curF = cp[facei] - 1;
//
//                    // Is the face on the boundary?
//                    if (curF >= mesh_.nInternalFaces())
//                    {
//                        label curBPatch =
//                            mesh_.boundaryMesh().whichPatch(curF);
//
//                        if (!patchFields(curBPatch))
//                        {
//                            patchFields.set
//                            (
//                                curBPatch,
//                                unallocatedFvPatchField<Type>::New
//                                (
//                                    mesh_.boundary()[curBPatch].type(),
//                                    mesh_.boundary()[curBPatch],
//                                    DimensionedField
//                                    <
//                                        Type,
//                                        unallocatedFvMesh
//                                    >::null()
//                                )
//                            );
//                        }
//
//                        // add the face
//                        label curPatchFace =
//                            mesh_.boundaryMesh()
//                                [curBPatch].whichFace(curF);
//
//                        patchFields[curBPatch][curPatchFace] =
//                            curProcPatch[facei];
//                    }
//                }
//            }
//        }
//    }
//
//    forAll(mesh_.boundary(), patchi)
//    {
//        // add empty patches
//        if
//        (
//            isType<emptyFvPatch>(mesh_.boundary()[patchi])
//         && !patchFields(patchi)
//        )
//        {
//            patchFields.set
//            (
//                patchi,
//                unallocatedFvPatchField<Type>::New
//                (
//                    emptyFvPatchField<Type>::typeName,
//                    mesh_.boundary()[patchi],
//                    DimensionedField<Type, unallocatedFvMesh>::null()
//                )
//            );
//        }
//    }
//
//
//    // Now construct and write the field
//    // setting the internalField and patchFields
//    return tmp
//    <GeometricField<Type, unallocatedFvPatchField, unallocatedFvMesh>>
//    (
//        new GeometricField<Type, unallocatedFvPatchField, unallocatedFvMesh>
//        (
//            fieldIoObject,
//            mesh_,
//            procFields[0].dimensions(),
//            internalField,
//            patchFields
//        )
//    );
//}


// ************************************************************************* //
