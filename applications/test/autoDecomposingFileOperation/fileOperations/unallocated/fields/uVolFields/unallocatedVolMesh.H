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
    Foam::unallocatedVolMesh

Description
    Adaptor to construct 'volume' GeometricFields on unallocatedFvMesh.

\*---------------------------------------------------------------------------*/

#ifndef unallocatedVolMesh_H
#define unallocatedVolMesh_H

#include "unallocatedFvMesh.H"
#include "GeoMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalMeshData;

/*---------------------------------------------------------------------------*\
                     Class unallocatedVolMesh Declaration
\*---------------------------------------------------------------------------*/

class unallocatedVolMesh
:
    public GeoMesh<unallocatedFvMesh>
{
public:

    // Constructors

        //- Construct from unallocatedFvMesh
        explicit unallocatedVolMesh(const unallocatedFvMesh& mesh)
        :
            GeoMesh<unallocatedFvMesh>(mesh)
        {}


    // Member Functions

        // volMesh functionality

            //- Return number of cells
            static label size(const Mesh& mesh)
            {
                return mesh.nCells();
            }

            //- Return reference to boundary mesh
            const unallocatedFvBoundaryMesh& boundary() const
            {
                return mesh_.boundary();
            }

            //- Return parallel info
            const globalMeshData& globalData() const
            {
                NotImplemented;
                return mesh_.globalData();
            }

            //- Return database
            const objectRegistry& thisDb() const
            {
                return mesh_.thisDb();
            }

            //- Return the top-level database
            const Time& time() const
            {
                return thisDb().time();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
