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
    Foam::unallocatedFvsMesh

Description
    Type of fvMesh without geometry storage.

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvsMesh_H
#define unallocatedFvsMesh_H

#include "unallocatedFvMesh.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class globalMeshData;

/*---------------------------------------------------------------------------*\
                           Class unallocatedFvsMesh Declaration
\*---------------------------------------------------------------------------*/

class unallocatedFvsMesh
:
    public unallocatedFvMesh,
    public surfaceMesh
{
////    // Permanent data
////
////        const word name_;
////
////        //- Database
////        const objectRegistry& db_;
////
////        //- Number of internal faces
////        const label nInternalFaces_;
////
////        //- Number of faces
////        label nFaces_;
////
////        //- Number of cells
////        const label nCells_;
////
////        //- Boundary mesh
////        unallocatedFvBoundaryMesh boundary_;
////
////        const globalMeshData& globalData_;
//
//
//    // Private Member Functions
//
//        //- Disallow default bitwise copy construct
//        unallocatedFvsMesh(const unallocatedFvsMesh&);
//
//        //- Disallow default bitwise assignment
//        void operator=(const unallocatedFvsMesh&);


public:

    // Declare name of the class and its debug switch
    ClassName("unallocatedFvsMesh");


    typedef unallocatedFvsMesh Mesh;
//    typedef unallocatedFvBoundaryMesh BoundaryMesh;


    // Constructors

        //- Construct from parts. nFaces, nInternalFaces derived from patch
        //  info.
        unallocatedFvsMesh
        (
            const word& name,
            const objectRegistry& db,
            const label nCells,
            const wordList& patchNames,
            const labelList& patchSizes,
            const labelList& patchStarts,
            const List<wordList>& patchGroups,
            const globalMeshData& globalData
        );


    //- Destructor
    ~unallocatedFvsMesh();


    // Member Functions

        // volMesh functionality

            //- Return number of cells
            static label size(const Mesh& mesh)
            {
                return mesh.nInternalFaces();
            }

//            //- Return reference to boundary mesh
//            const unallocatedFvBoundaryMesh& boundary() const
//            {
//                return boundary_;
//            }
//
//            //- Return parallel info
//            const globalMeshData& globalData() const
//            {
//                NotImplemented;
//                return globalData_;
//            }
//
//            //- Return the top-level database
//            const Time& time() const
//            {
//                return db_.time();
//            }
//
//            //- Return database
//            const objectRegistry& thisDb() const
//            {
//                return db_;
//            }


//        // Basic mesh info
//
//            //- Override database name since this might be Time
//            const word& name() const
//            {
//                return name_;
//            }
//
//            label nCells() const
//            {
//                return nCells_;
//            }
//
//            label nFaces() const
//            {
//                return nFaces_;
//            }
//
//            label nInternalFaces() const
//            {
//                return nInternalFaces_;
//            }


//    // Member Operators
//
//        bool operator!=(const unallocatedFvsMesh&) const;
//        bool operator==(const unallocatedFvsMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
