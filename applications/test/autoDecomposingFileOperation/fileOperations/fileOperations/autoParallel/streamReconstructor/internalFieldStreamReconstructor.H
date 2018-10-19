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
    Foam::internalFieldStreamReconstructor

Description

SourceFiles
    internalFieldStreamReconstructor.C
    internalFieldStreamReconstructors.C

\*---------------------------------------------------------------------------*/

#ifndef internalFieldStreamReconstructor_H
#define internalFieldStreamReconstructor_H

#include "streamReconstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class internalFieldStreamReconstructor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class internalFieldStreamReconstructor
:
    public streamReconstructor
{
public:

    // Constructors

        //- Construct from components
        internalFieldStreamReconstructor()
        {}


    //- Destructor
    virtual ~internalFieldStreamReconstructor()
    {}


    // Member Functions

        //- Reconstruct and stream
        virtual bool reconstruct
        (
            const IOobject& io,
            const bool isOriented,      // apply flip (if available)
            Ostream& os
        ) const;

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
        ) const
        {
            NotImplemented;
            return false;
        }

        //- Parallel reconstruct and (on master) write
        virtual bool reconstruct
        (
            const parUnallocatedFvFieldReconstructor& reconstructor,
            const regIOobject& thisIO,
            const bool isOriented,      // apply flip (if available)
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const
        {
            NotImplemented;
            return false;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "internalFieldStreamReconstructor.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
