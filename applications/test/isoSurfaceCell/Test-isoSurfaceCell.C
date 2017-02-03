/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

Application
    Test-isoSurfaceCell

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "fvMesh.H"
#include "isoSurfaceCell.H"
#include "OBJstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"


    scalarField ccx(mesh.cellCentres().component(vector::Y));
    scalarField pointsx(mesh.points().component(vector::Y));

    const scalar value = 1.1*average(pointsx);

    isoSurfaceCell iso
    (
        mesh,
        ccx,
        pointsx,
        value,
        false       // regularise,
    );

    Pout<< "iso:" << iso << endl;
    Pout<< "usesCellCentre:" << iso.usesCellCentre() << endl;

    iso.write("iso.obj");

    OBJstream str("ccPoints.obj");
    forAll(iso.usesCellCentre(), pointi)
    {
        if (iso.usesCellCentre()[pointi])
        {
            str.write(iso.points()[pointi]);
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
