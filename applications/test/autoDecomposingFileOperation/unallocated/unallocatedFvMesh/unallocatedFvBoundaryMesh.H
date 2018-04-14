/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::unallocatedFvBoundaryMesh

Description

SourceFiles
    //unallocatedFvBoundaryMeshI.H
    //unallocatedFvBoundaryMesh.C
    //unallocatedFvBoundaryMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvBoundaryMesh_H
#define unallocatedFvBoundaryMesh_H

#include "PtrList.H"
#include "HashSet.H"
#include "unallocatedGenericFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class unallocatedFvBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class unallocatedFvBoundaryMesh
:
    public PtrList<unallocatedGenericFvPatch>
{
    // Private data

        mutable autoPtr<HashTable<labelList, word>> groupPatchIDsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        unallocatedFvBoundaryMesh(const unallocatedFvBoundaryMesh&);

        //- Disallow default bitwise assignment
        void operator=(const unallocatedFvBoundaryMesh&);


public:

    // Constructors

        //- Construct null
        unallocatedFvBoundaryMesh()
        :
            PtrList<unallocatedGenericFvPatch>()
        {}


    // Member functions

        //- Find patch index given a name
        label findPatchID(const word& patchName) const;

        //- Find patch indices given a name
        labelList findIndices(const keyType&, const bool useGroups) const;

        //- Return patch index for a given face label
        label whichPatch(const label faceIndex) const;

        //- Per patch group the patch indices
        const HashTable<labelList, word>& groupPatchIDs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
