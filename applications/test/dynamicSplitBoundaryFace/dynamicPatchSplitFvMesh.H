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
    Foam::dynamicPatchSplitFvMesh

Description
    A fvMesh with built-in refinement.

SourceFiles
    dynamicPatchSplitFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicPatchSplitFvMesh_H
#define dynamicPatchSplitFvMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicPatchSplitFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicPatchSplitFvMesh
:
    public dynamicFvMesh
{
private:

        //- Disallow default bitwise copy construct
        dynamicPatchSplitFvMesh(const dynamicPatchSplitFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicPatchSplitFvMesh&);

public:

    //- Runtime type information
    TypeName("dynamicPatchSplitFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicPatchSplitFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicPatchSplitFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
