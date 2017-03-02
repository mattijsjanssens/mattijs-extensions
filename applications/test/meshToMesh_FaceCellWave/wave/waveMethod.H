/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::waveMethod

Description
    Direct (one-to-one cell correspondence) mesh-to-mesh interpolation class

SourceFiles
    waveMethod.C

\*---------------------------------------------------------------------------*/

#ifndef waveMethod_H
#define waveMethod_H

#include "meshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class waveMethod Declaration
\*---------------------------------------------------------------------------*/

class waveMethod
:
    public meshToMeshMethod
{
    // Private member functions

        //- Calculate addressing and weights
        void calculate
        (
            const polyMesh& src,
            const polyMesh& tgt,
            labelListList& srcToTgtAddr,
            scalarListList& srcToTgtWght
        ) const;

        //- Disallow default bitwise copy construct
        waveMethod(const waveMethod&);

        //- Disallow default bitwise assignment
        void operator=(const waveMethod&);


public:

    //- Run-time type information
    TypeName("wave");

    //- Construct from source and target meshes
    waveMethod(const polyMesh& src, const polyMesh& tgt);

    //- Destructor
    virtual ~waveMethod();


    // Member Functions

        // Evaluate

            //- Calculate addressing and weights
            virtual void calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                labelListList& tgtToTgtAddr,
                scalarListList& tgtToTgtWght
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
