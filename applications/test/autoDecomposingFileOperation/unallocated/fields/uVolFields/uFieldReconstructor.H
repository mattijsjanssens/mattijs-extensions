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
    Foam::uFieldReconstructor

Description
    MeshObject wrapper around data for reconstructing uVolFields.

SourceFiles
    uFieldReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef uFieldReconstructor_H
#define uFieldReconstructor_H

//#include "MeshObject.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class unallocatedFvMesh;
class unallocatedFvFieldReconstructor;


/*---------------------------------------------------------------------------*\
                     Class uFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class uFieldReconstructor
// :
//     public MeshObject<polyMesh, TopologicalMeshObject, uFieldReconstructor>
:
    public regIOobject
{
    // Private Member Data

        mutable autoPtr<unallocatedFvMesh> baseMeshPtr_;
        mutable PtrList<Time> procDatabases_;
        mutable PtrList<labelIOList> cellProcAddressing_;
        mutable PtrList<labelIOList> faceProcAddressing_;
        mutable PtrList<labelIOList> boundaryProcAddressing_;
        mutable PtrList<unallocatedFvMesh> procMeshes_;

        mutable autoPtr<unallocatedFvFieldReconstructor> reconstructorPtr_;


    // Private Member Functions

        //- Given parent IOobject read the processor databases
        void readProcDatabases(const IOobject& io, const label nProcs) const;

        //- Read processor meshes and addressing
        void readProcMeshes(const Time&, const fileName& instance) const;


public:

    // Declare name of the class and its debug switch
    TypeName("uFieldReconstructor");


    // Constructors

        //- Constructor given objectRegistry (undecomposed mesh). Registry
        //  itself not used apart from its path and time
        explicit uFieldReconstructor(const objectRegistry& obr);

        //- Lookup or create
        static const uFieldReconstructor& New(const objectRegistry& obr);


    //- Destructor
    virtual ~uFieldReconstructor()
    {}


    // Member Functions

        //- Undecomposed mesh
        const unallocatedFvMesh& baseMesh() const
        {
            return baseMeshPtr_();
        }

        //- Decomposed meshes
        const PtrList<unallocatedFvMesh>& procMeshes() const
        {
            return procMeshes_;
        }

        //- Field reconstruction
        const unallocatedFvFieldReconstructor& reconstructor() const
        {
            return reconstructorPtr_();
        }

        virtual bool writeData(Foam::Ostream&) const
        {
            NotImplemented;
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
