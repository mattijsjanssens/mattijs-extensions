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
    Foam::lagrangianStreamReconstructor

Description

SourceFiles
    lagrangianStreamReconstructor.C
    lagrangianStreamReconstructorNew.C

\*---------------------------------------------------------------------------*/

#ifndef lagrangianStreamReconstructor_H
#define lagrangianStreamReconstructor_H

#include "streamReconstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class lagrangianStreamReconstructor Declaration
\*---------------------------------------------------------------------------*/

class lagrangianStreamReconstructor
:
    public streamReconstructor
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        lagrangianStreamReconstructor(const lagrangianStreamReconstructor&);

        //- Disallow default bitwise assignment
        void operator=(const lagrangianStreamReconstructor&);


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            lagrangianStreamReconstructor,
            cloudName,
            (
                const word& cloudName
            ),
            (cloudName)
        );


    // Constructors

        //- Construct null
        lagrangianStreamReconstructor()
        {}


    // Selectors

        //- Return a pointer to a new cloud-field reconstructor created on
        //  freestore from components
        static autoPtr<lagrangianStreamReconstructor> New
        (
            const word& type,
            const word& cloudName
        );


    //- Destructor
    virtual ~lagrangianStreamReconstructor()
    {}


    // Member Functions

        //- (on master) Read all processor's field, reconstruct and stream
        virtual bool reconstruct
        (
            const IOobject& io,
            const bool isOriented,      // apply flip (if available)
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
            const bool isOriented,      // apply flip (if available)
            Ostream& os
        ) const = 0;

        //- Parallel reconstruct and (on master) write
        virtual bool reconstruct
        (
            const parUnallocatedFvFieldReconstructor& reconstructor,
            const regIOobject& thisIO,
            const bool isOriented,      // apply flip (if available)
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
