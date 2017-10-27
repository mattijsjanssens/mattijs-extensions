/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::unallocatedFvMesh

Description
    Mesh representing a set of points created from polyMesh.

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvMesh_H
#define unallocatedFvMesh_H

//#include "GeoMesh.H"
//#include "MeshObject.H"
//#include "polyMesh.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class unallocatedFvBoundaryMesh;
class globalMeshData;

/*---------------------------------------------------------------------------*\
                           Class unallocatedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class unallocatedFvMesh
:
    public volMesh
{
    // Permanent data

        //- Database
        const objectRegistry& db_;

        //- Size of internal
        const label nCells_;

        //- Boundary mesh
        const unallocatedFvBoundaryMesh& boundary_;

        const globalMeshData& globalData_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        unallocatedFvMesh(const unallocatedFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const unallocatedFvMesh&);


public:

    // Declare name of the class and its debug switch
    ClassName("unallocatedFvMesh");


    typedef unallocatedFvMesh Mesh;
    typedef unallocatedFvBoundaryMesh BoundaryMesh;


    // Constructors

        //- Construct from polyMesh
        explicit unallocatedFvMesh
        (
            //const fvMesh& procMesh,
            const objectRegistry& db,
            const label nCells,
            const unallocatedFvBoundaryMesh& boundary,
            const globalMeshData& globalData
        );


    //- Destructor
    ~unallocatedFvMesh();


    // Member Functions

        //- Return number of cells
        static label size(const Mesh& mesh)
        {
            return mesh.nCells_;
        }

        //- Return reference to boundary mesh
        const unallocatedFvBoundaryMesh& boundary() const
        {
            return boundary_;
        }

        //- Return parallel info
        const globalMeshData& globalData() const
        {
            return globalData_;
        }

        //- Return the top-level database
        const Time& time() const
        {
            return db_.time();
        }

        //- Return database. For now is its polyMesh.
        const objectRegistry& thisDb() const
        {
            return db_;
        }


    // Member Operators

        bool operator!=(const unallocatedFvMesh&) const;
        bool operator==(const unallocatedFvMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
