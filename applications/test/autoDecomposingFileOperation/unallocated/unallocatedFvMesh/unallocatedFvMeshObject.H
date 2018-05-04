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
    Foam::unallocatedFvMeshObject

Description
    Construct unallocatedFvMesh from fvMesh.

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvMeshObject_H
#define unallocatedFvMeshObject_H

#include "fvMesh.H"
#include "unallocatedFvMesh.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class unallocatedFvMeshObject Declaration
\*---------------------------------------------------------------------------*/

class unallocatedFvMeshObject
:
    public MeshObject<fvMesh, TopologicalMeshObject, unallocatedFvMeshObject>,
    public unallocatedFvMesh
{
    // Private Member Functions

        static wordList patchNames(const fvMesh&);
        static labelList patchSizes(const fvMesh&);
        static labelList patchStarts(const fvMesh&);
        static List<wordList> patchGroups(const fvMesh&);


public:

    TypeName("unallocatedFvMesh");

    // Constructors

        //- Construct from fvMesh
        unallocatedFvMeshObject(const fvMesh& mesh);


    //- Destructor
    virtual ~unallocatedFvMeshObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
