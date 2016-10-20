/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::ACMIWallFvPatch

Description
    Foam::ACMIWallFvPatch

SourceFiles
    ACMIWallFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef ACMIWallFvPatch_H
#define ACMIWallFvPatch_H

#include "wallFvPatch.H"
#include "ACMIWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ACMIWallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class ACMIWallFvPatch
:
    public wallFvPatch
{

public:

    //- Runtime type information
    TypeName(ACMIWallPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        ACMIWallFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            wallFvPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
