/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::fileFormats::ensightMeshReader

Description

Notes

SourceFiles
    ensightMeshReader.C

\*---------------------------------------------------------------------------*/

#ifndef ensightMeshReader_H
#define ensightMeshReader_H

#include "meshReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::ensightMeshReader Declaration
\*---------------------------------------------------------------------------*/

class ensightMeshReader
:
    public meshReader
{
    // Private data

        //- Retain solid cell types
        //bool keepSolids_;


protected:

    // Protected Data

        //- Point labels (imported Point numbering not necessarily contiguous)
        labelList mapToFoamPointId_;

        //- Cell labels (imported Cell numbering not necessarily contiguous)
        labelList mapToFoamCellId_;


    // Protected Member Functions

        //- Rotate face so lowest vertex is first
        const face& rotateFace
        (
            const face& f,
            face& rotatedFace
        ) const;

        //- Read the mesh from the file(s)
        virtual bool readGeometry(const scalar scaleFactor = 1.0);


public:

    //- Runtime type information
    TypeName("ensightMeshReader");


    // Constructors

        //- Construct from case name
        ensightMeshReader
        (
            const fileName& prefix,
            const objectRegistry& registry,
            const scalar scaleFactor = 1.0,
            const bool keepSolids = false
        );


    //- Destructor
    virtual ~ensightMeshReader() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
