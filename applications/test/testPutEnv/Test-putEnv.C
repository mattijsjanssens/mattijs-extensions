/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Mattijs Janssens
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

#include "fvMesh.H"
#include "volFields.H"
#include "argList.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

static char buf[200];


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    string s("VAR=");
    strncpy(buf, s.c_str(), 200-1);
    putenv(buf);

    for (label i=0; i < 1000000; i++)
    {
        string s("VAR=" + Foam::name(i));
        strncpy(buf, s.c_str(), 200-1);

        //Pout<< "buf:" << buf << endl;
        //string contents(getEnv("VAR"));
        //DebugVar(contents);
    }

    string contents(getEnv("VAR"));
    DebugVar(contents);

    return 0;
}


// ************************************************************************* //
