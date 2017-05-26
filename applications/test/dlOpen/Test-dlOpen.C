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
    Test-dlOpen

Description
    Test loading and unloading of libraries

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "argList.H"
//#include "Time.H"
#include "stringListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("libraries");
    #include "setRootCase.H"

    fileNameList fNames
    (
        readDir
        (
            getEnv("FOAM_LIBBIN"),
            fileName::FILE,
            false
        )
    );
    fNames = UIndirectList<fileName>(fNames, findStrings(".*\\.so", fNames))();
    //fNames.setSize(1, "libengine.so");

    //DebugVar(args.arg(1));
    //fileNameList fNames(IStringStream(args.arg(1))());

    Info<< "fNames:" << fNames << endl;

    forAll(fNames, i)
    {
        const fileName& fName = fNames[i];

        Info<< "Loading " << fName << endl;
        void* handle = dlOpen(fName);
        if (!handle)
        {
            WarningInFunction << "Cannot dlOpen " << fName << endl;
        }
        else
        {
            Info<< "Unloading " << fName << endl;
            bool ok = dlClose(handle);
            if (!ok)
            {
               WarningInFunction << "Cannot dlClose " << fName << endl;
            }
        }
    }

    Info<< "end" << endl;
}


// ************************************************************************* //
