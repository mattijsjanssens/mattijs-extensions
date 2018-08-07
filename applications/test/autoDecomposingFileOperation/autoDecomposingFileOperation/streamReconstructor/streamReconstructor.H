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
    Foam::streamReconstructor

Description

SourceFiles
    streamReconstructor.C
    streamReconstructorNew.C

\*---------------------------------------------------------------------------*/

#ifndef streamReconstructor_H
#define streamReconstructor_H

#include "fvMesh.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class unallocatedFvMesh;
class parUnallocatedFvFieldReconstructor;

/*---------------------------------------------------------------------------*\
                     Class streamReconstructor Declaration
\*---------------------------------------------------------------------------*/

class streamReconstructor
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        streamReconstructor(const streamReconstructor&);

        //- Disallow default bitwise assignment
        void operator=(const streamReconstructor&);


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            streamReconstructor,
            word,
            (
            ),
            ()
        );


    // Constructors

        //- Construct null
        streamReconstructor()
        {}


    // Selectors

        //- Return a pointer to a new patch created on freestore from
        //  components
        static autoPtr<streamReconstructor> New(const word& type);


    //- Destructor
    virtual ~streamReconstructor()
    {}


    // Member Functions

        //- (on master) Read all processor's field, reconstruct and stream
        virtual bool reconstruct
        (
            const fvMesh& mesh,
            const IOobject& io,
            Ostream& os
        ) const = 0;

        //- Read, (in parallel) decompose and stream
        virtual bool decompose
        (
            const parUnallocatedFvFieldReconstructor& reconstructor,
            const unallocatedFvMesh& baseMesh,
            const IOobject& baseIO,

            const unallocatedFvMesh& thisMesh,
            const IOobject& thisIO,
            Ostream& os
        ) const = 0;

        //- Parallel reconstruct and (on master) write
        virtual bool reconstruct
        (
            const parUnallocatedFvFieldReconstructor& reconstructor,
            const regIOobject& thisIO,
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
