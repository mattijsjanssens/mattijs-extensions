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
    Test-parallel-communicators

Description
    Checks communication using user-defined communicators

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "zstr.hpp"
#include "OFstream.H"
#include "OStringStream2.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    //#include "createTime.H"

    std::stringbuf* bufPtr = new std::stringbuf;

    // Stream data into buffer
    {
        zstr::ostream os(bufPtr);
        //ostream os(bufPtr);
        os << "BLABLA" << std::endl;
        cout<< "os:" << bufPtr->str() << std::endl;
    }

    //// Write buffer into OFstream
    //{
    //    OFstream os("file_working.gz", IOstream::streamFormat::BINARY);
    //    os.stdStream() << bufPtr->str();
    //}
    {
        OFstream os("file.gz", IOstream::streamFormat::BINARY);
        os.writeQuoted(bufPtr->str(), false);
    }


    // Retrieve data from buffer
    {
        zstr::istream is(bufPtr);
        std::string a;
        is >> a;
        cout<< "a:" << a << std::endl;
    }

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
