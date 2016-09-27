/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Test-List

Description
    Simple tests and examples of use of List

See also
    Foam::List

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "Map.H"
#include "SortableList.H"
#include "Field.H"

#include "IOstreams.H"
#include "scalar.H"
#include "vector.H"
#include "mergePoints.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Map<label> m = {{200, 1},{300, 2}};
    SortableList<label> l({300, 100, 200});
    DebugVar(l);

    scalarField a(3);
    a[0] = 0.0;
    a[1] = 1.1;
    a[2] = 2.2;

    labelList addressing(2);
    addressing[0] = 2;
    addressing[1] = 1;

    UIndirectList<scalar> b(a, addressing);
    DebugVar(b);

//    Field<scalar> c(b);
//    DebugVar(c);

    return 0;
}

// ************************************************************************* //
