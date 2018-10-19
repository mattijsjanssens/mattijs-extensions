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
    Foam::FieldStreamReconstructor

Description

SourceFiles
    FieldStreamReconstructor.C
    FieldStreamReconstructors.C

\*---------------------------------------------------------------------------*/

#ifndef FieldStreamReconstructor_H
#define FieldStreamReconstructor_H

#include "lagrangianStreamReconstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class FieldStreamReconstructor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FieldStreamReconstructor
:
    public lagrangianStreamReconstructor
{
public:

    // Declare name of the class (used in the run-time selection)
    ClassNameNoDebug("scalarField");

    // Constructors

        //- Construct from components
        FieldStreamReconstructor(const word& cloudName)
        :
            lagrangianStreamReconstructor()
        {}


    //- Destructor
    virtual ~FieldStreamReconstructor()
    {}


    // Member Functions

        //- Reconstruct and stream
        virtual bool reconstruct
        (
            //const fvMesh& mesh,
            const IOobject& io,
            const bool isOriented,
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
            const bool isOriented,
            Ostream& os
        ) const;

        virtual bool reconstruct
        (
            const parUnallocatedFvFieldReconstructor& reconstructor,
            const regIOobject& thisIO,
            const bool isOriented,
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FieldStreamReconstructor.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
