/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
#include "IOobjectList.H"
#include "mapDistributePolyMesh.H"
#include "processorFvPatch.H"

#include "directFvPatchFieldMapper.H"
#include "distributedUnallocatedDirectFieldMapper.H"
#include "distributedUnallocatedDirectFvPatchFieldMapper.H"

#include "unallocatedFvMesh.H"
#include "unallocatedFvBoundaryMesh.H"

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


// Reconstruct a field onto the baseMesh
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
Foam::parUnallocatedFvFieldReconstructor::reconstructFvVolumeField
(
    const GeometricField<Type, unallocatedFvPatchField, unallocatedVolMesh>& fld
) const
{
    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap()
    );

    Field<Type> internalField(fld.internalField(), mapper);



    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on mesh, not baseMesh

    PtrList<unallocatedFvPatchField<Type> > patchFields
    (
        fld.mesh().boundary().size()
    );

    const typename GeometricField
    <
        Type,
        unallocatedFvPatchField,
        unallocatedVolMesh
    >::Boundary& bfld = fld.boundaryField();

    forAll(bfld, patchI)
    {
        if (patchFaceMaps_.set(patchI))
        {
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


    PtrList<unallocatedFvPatchField<Type> > basePatchFields
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

            const unallocatedFvPatchField<Type>& pfld = patchFields[patchI];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            basePatchFields.set
            (
                patchI,
                unallocatedFvPatchField<Type>::New
                (
                    pfld,
                    basePatch,
                    DimensionedField<Type, unallocatedVolMesh>::null(),
                    dummyMapper
                )
            );
        }
    }

    // Add some empty patches on remaining patches (tbd.probably processor
    // patches)
    forAll(basePatchFields, patchI)
    {
        if (patchI >= patchFields.size() || !patchFields.set(patchI))
        {
            basePatchFields.set
            (
                patchI,
                unallocatedFvPatchField<Type>::New
                (
                    emptyFvPatchField<Type>::typeName,
                    baseMesh_.boundary()[patchI],
                    DimensionedField<Type, unallocatedVolMesh>::null()
                )
            );
        }
    }

    // Construct a volField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_.thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    return tmp
    <
        GeometricField
        <
            Type,
            unallocatedFvPatchField,
            unallocatedVolMesh
        >
    >
    (
        new GeometricField<Type, unallocatedFvPatchField, unallocatedVolMesh>
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );
}


// ************************************************************************* //
