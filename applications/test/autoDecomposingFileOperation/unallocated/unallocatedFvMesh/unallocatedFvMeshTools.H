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
    Foam::unallocatedFvMeshTools

Description
    A collection of tools for operating on an unallocatedFvMesh.

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvMeshTools_H
#define unallocatedFvMeshTools_H

#include "unallocatedFvMesh.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class unallocatedFvMeshTools Declaration
\*---------------------------------------------------------------------------*/

class unallocatedFvMeshTools
{
    //- Send local map to master and convert into mapDistribute
    static autoPtr<mapDistribute> calcReconstructMap
    (
        const labelList& localToGlobal,
        const bool localHasFlip
    );


public:

    //- Read the *procAddressing files (from provided location) and
    //  construct a single mesh reconstruction map from it
    static autoPtr<mapDistributePolyMesh> readReconstructMap
    (
        const IOobject& io
    );

    //- Read boundary (searches for "faces" from provided instance) and
    //  construct placeholder mesh. Has illegal references to underlying
    //  globalMeshData (hopefully unused)
    static autoPtr<unallocatedFvMesh> newMesh
    (
        const IOobject& io,
        const label nCells
    );

    //- Read owner&boundary (searches for "faces" from provided instance) and
    //  construct placeholder mesh. Has illegal references to underlying
    //  globalMeshData (hopefully unused)
    static autoPtr<unallocatedFvMesh> newMesh(const IOobject& io);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
