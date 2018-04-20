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
    Foam::unallocatedFvMesh

Description
    Type of fvMesh without geometry storage.

\*---------------------------------------------------------------------------*/

#ifndef unallocatedFvMesh_H
#define unallocatedFvMesh_H

#include "unallocatedFvBoundaryMesh.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalMeshData;
//class objectRegistry;

/*---------------------------------------------------------------------------*\
                      Class unallocatedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class unallocatedFvMesh
{
    // Private Member Data

        const word name_;

        //- Database
        const objectRegistry& db_;

        //- Number of internal faces
        const label nInternalFaces_;

        //- Number of faces
        label nFaces_;

        //- Number of cells
        const label nCells_;

        //- Boundary mesh
        unallocatedFvBoundaryMesh boundary_;

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

        //- Construct from parts. nFaces, nInternalFaces derived from patch
        //  info.
        unallocatedFvMesh
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
    ~unallocatedFvMesh();


    // Member Functions

        // Basic mesh info

            //- Return parallel info
            const globalMeshData& globalData() const
            {
                NotImplemented;
                return globalData_;
            }

            //- Return database
            const objectRegistry& thisDb() const
            {
                return db_;
            }

            //- Return the top-level database
            const Time& time() const
            {
                return db_.time();
            }

            //- Override database name since this might be Time
            const word& name() const
            {
                return name_;
            }

            label nCells() const
            {
                return nCells_;
            }

            label nFaces() const
            {
                return nFaces_;
            }

            label nInternalFaces() const
            {
                return nInternalFaces_;
            }

            //- Return reference to boundary mesh
            const unallocatedFvBoundaryMesh& boundary() const
            {
                return boundary_;
            }


        //- Return info proxy.
        //  Used to print mesh information to a stream
        InfoProxy<unallocatedFvMesh> info() const
        {
            return *this;
        }


    // Member Operators

        bool operator!=(const unallocatedFvMesh&) const;
        bool operator==(const unallocatedFvMesh&) const;
};


template<>
Ostream& operator<<(Ostream&, const InfoProxy<unallocatedFvMesh>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
