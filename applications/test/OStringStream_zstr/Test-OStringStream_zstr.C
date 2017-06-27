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
#include "IStringStream2.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    //#include "createTime.H"

    {
        autoPtr<std::stringbuf> bufPtr(new std::stringbuf);

        // Stream data into buffer
        {
            zstr::ostream os(bufPtr.operator->());
            //ostream os(bufPtr);
            os << "BLABLA" << std::endl;
            cout<< "os:" << bufPtr->str() << std::endl;
        }

        // Write buffer into OFstream
        {
            OFstream os("file.gz", IOstream::streamFormat::BINARY);
            //os.stdStream() << bufPtr->str();
            os.writeQuoted(bufPtr->str(), false);
        }

        // Retrieve data from buffer
        {
            zstr::istream is(bufPtr.operator->());
            std::string a;
            is >> a;
            cout<< "a:" << a << std::endl;
        }
    }

    // Wrapped up in OStringStream
    string testStr;
    {
        // Test compressed stream
        OStringStream2 os
        (
            IOstream::streamFormat::ASCII,
            IOstream::currentVersion,
            IOstream::compressionType::UNCOMPRESSED
        );
        os << "MORE" << endl;
        os << "BLA" << endl;

        // Test copy
        OStringStream2 os2(os);

        testStr = os2.str();

        DebugVar(testStr);
    }

    // Read from compressed string
    {
        IStringStream2 is
        (
            testStr.c_str(),        //testStr,
            IOstream::streamFormat::ASCII,
            IOstream::currentVersion,
            IOstream::compressionType::UNCOMPRESSED
        );
        word word1(is);
        DebugVar(word1);
        word word2(is);
        DebugVar(word2);
    }


    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
