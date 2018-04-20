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

Class
    Foam::unallocatedFvFieldReconstructor

Description
    Finite volume reconstructor for (unallocated) volume fields.

SourceFiles
    unallocatedFvFieldReconstructor.C
    unallocatedFvFieldReconstructorReconstructFields.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvFieldReconstructor_H
#define unallocatedFvFieldReconstructor_H

#include "unallocatedVolMesh.H"
#include "unallocatedFvPatchField.H"
#include "labelIOList.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class unallocatedFvFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class unallocatedFvFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        unallocatedFvMesh& mesh_;

        //- List of processor meshes
        const PtrList<unallocatedFvMesh>& procMeshes_;

        //- List of processor face addressing lists
        const PtrList<labelIOList>& faceProcAddressing_;

        //- List of processor cell addressing lists
        const PtrList<labelIOList>& cellProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;

        //- Number of fields reconstructed
        label nReconstructed_;


    // Private Member Functions

        //- Find entry with same key
        label findEntry
        (
            const List<Pair<word>>& lst,
            const word& key
        ) const;


        //- Disallow default bitwise copy construct
        unallocatedFvFieldReconstructor
        (
            const unallocatedFvFieldReconstructor&
        );

        //- Disallow default bitwise assignment
        void operator=(const unallocatedFvFieldReconstructor&);


public:

    // Constructors

        //- Construct from components
        unallocatedFvFieldReconstructor
        (
            unallocatedFvMesh& mesh,
            const PtrList<unallocatedFvMesh>& procMeshes,
            const PtrList<labelIOList>& faceProcAddressing,
            const PtrList<labelIOList>& cellProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );

    // Member Functions

        //- Fix generic non-uniform reading
        template<class Type>
        void fixGenericNonuniform
        (
            PtrList
            <
                GeometricField
                <
                    Type,
                    unallocatedFvPatchField,
                    unallocatedVolMesh
                >
            >& procFields
        ) const;

        //- Reconstruct volume internal field
        template<class Type>
        tmp<DimensionedField<Type, unallocatedVolMesh>>
        reconstructFvVolumeInternalField
        (
            const IOobject& fieldIoObject,
            const PtrList<DimensionedField<Type, unallocatedVolMesh>>& pFields
        ) const;

        //- Reconstruct volume field
        template<class Type>
        tmp<GeometricField<Type, unallocatedFvPatchField, unallocatedVolMesh>>
        reconstructFvVolumeField
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
            >&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "unallocatedFvFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
