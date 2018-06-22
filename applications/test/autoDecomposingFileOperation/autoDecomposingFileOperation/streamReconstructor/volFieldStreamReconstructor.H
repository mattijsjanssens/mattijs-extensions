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
    Foam::volFieldStreamReconstructor

Description

SourceFiles
    volFieldStreamReconstructor.C
    volFieldStreamReconstructors.C

\*---------------------------------------------------------------------------*/

#ifndef volFieldStreamReconstructor_H
#define volFieldStreamReconstructor_H

#include "streamReconstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class volFieldStreamReconstructor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class volFieldStreamReconstructor
:
    public streamReconstructor
{
public:

    // Constructors

        //- Construct from components
        volFieldStreamReconstructor()
        {}


    //- Destructor
    virtual ~volFieldStreamReconstructor()
    {}


    // Member Functions

        //- Reconstruct and stream
        virtual bool reconstruct
        (
            const fvMesh& mesh,
            const IOobject& io,
            Ostream& os
        ) const;

        //- Read, decompose and stream
        virtual bool decompose
        (
            const parUnallocatedFvFieldReconstructor& reconstructor,
            const unallocatedFvMesh& baseMesh,
            const IOobject& baseIO,

            const unallocatedFvMesh& thisMesh,
            const IOobject& thisIO,
            Ostream& os
        ) const;

        virtual bool reconstruct
        (
            const parUnallocatedFvFieldReconstructor& reconstructor,
            const regIOobject& thisIO,
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volFieldStreamReconstructor.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
