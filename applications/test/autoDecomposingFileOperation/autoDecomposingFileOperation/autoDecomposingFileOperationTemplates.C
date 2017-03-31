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

\*---------------------------------------------------------------------------*/

#include "fvFieldDecomposer.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoField>
bool Foam::fileOperations::autoDecomposingFileOperation::decomposeAndWrite
(
    const IOobject& procIO,
    const IOobject& parentIO,
    const word& type,
    Ostream& os
) const
{
    if (type == GeoField::typeName)
    {
        Info<< "Decomposing field " << parentIO.name() << endl;

        const fvMesh& undecomposedMesh = baseMesh(parentIO.time());
        GeoField parentFld(parentIO, undecomposedMesh);

        const fvFieldDecomposer& volDecomposer(decomposer(procIO));

        // Decompose field and transfer to stream
        tmp<GeoField> fld(volDecomposer.decomposeField(parentFld));
        fld().writeData(os);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
