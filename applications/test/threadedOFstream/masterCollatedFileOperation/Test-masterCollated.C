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

Application
    Test-masterCollated

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "OFstreamCollator.H"
#include "masterUncollatedFileOperation.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"


    fileOperations::masterUncollatedFileOperation::maxMasterFileBufferSize =
        1000000000;
    const off_t maxBufferSize(100);
    const label size = 10;

    OFstreamCollator writer(maxBufferSize);

    string data(size, ' ');
    forAll(data, i)
    {
       data[i] = 'A' + Pstream::myProcNo();
    }

    bool ok = writer.write
    (
        string::typeName,
        runTime.path()/"myFile.txt",
        data,

        IOstream::BINARY,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED,
        false                       // append
    );


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
