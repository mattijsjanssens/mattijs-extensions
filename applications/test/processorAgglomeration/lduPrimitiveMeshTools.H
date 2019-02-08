/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::lduPrimitiveMeshTools

Description
    A collection of tools on lduMeshes.

SourceFiles
    lduPrimitiveMeshTools.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMeshTools_H
#define lduPrimitiveMeshTools_H

#include "lduPrimitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class lduPrimitiveMeshTools Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMeshTools
{
public:

    //- Write mesh (move to lduMesh? Needs IO functionality on lduInterface!)
    static bool writeData(const lduMesh& mesh, Ostream& os);

    //- Read mesh (move to lduPrimitiveMesh?)
    static autoPtr<lduPrimitiveMesh> readData(Istream& is);

    //- Exchange label data on interfaces
    static void swapCellData
    (
        const lduInterfacePtrsList& ifs,
        const labelUList& cellData,
        PtrList<labelField>& nbrData
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
