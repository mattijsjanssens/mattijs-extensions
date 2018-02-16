/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::parFvFieldReconstructor

Description
    Finite volume reconstructor for volume and surface fields.

    Runs in parallel. Reconstructs from procMesh to baseMesh. baseMesh
    is non-zero cells on processor0 only.

SourceFiles
    parFvFieldReconstructor.C
    parFvFieldReconstructorReconstructFields.C

\*---------------------------------------------------------------------------*/

#ifndef parFvFieldReconstructor_H
#define parFvFieldReconstructor_H

#include "PtrList.H"
#include "fvMesh.H"
#include "unallocatedFvMesh.H"
#include "unallocatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistributePolyMesh;
class mapDistributeBase;

/*---------------------------------------------------------------------------*\
                   Class parFvFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class parFvFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        unallocatedFvMesh& baseMesh_;

        //- Processor mesh reference
        const fvMesh& procMesh_;

        //- Distribution map reference
        const mapDistributePolyMesh& distMap_;

        //- Do I need to to write (usually master only)
        const bool isWriteProc_;

        //- Patch mappers
        PtrList<mapDistributeBase> patchFaceMaps_;


    // Private Member Functions

        //- Construct per-patch addressing
        void createPatchFaceMaps();

        //- Disallow default bitwise copy construct
        parFvFieldReconstructor(const parFvFieldReconstructor&);

        //- Disallow default bitwise assignment
        void operator=(const parFvFieldReconstructor&);


public:

    // Constructors

        //- Construct from components
        parFvFieldReconstructor
        (
            unallocatedFvMesh& baseMesh,
            const fvMesh& procMesh,
            const mapDistributePolyMesh& distMap,
            const bool isWriteProc
        );


    // Member Functions

        //- Reconstruct volume internal field
        template<class Type>
        tmp<DimensionedField<Type, unallocatedFvMesh> >
        reconstructFvVolumeInternalField
        (
            const DimensionedField<Type, volMesh>&
        ) const;

        //- Reconstruct volume field
        template<class Type>
        tmp<GeometricField<Type, unallocatedFvPatchField, unallocatedFvMesh> >
        reconstructFvVolumeField
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "parFvFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
