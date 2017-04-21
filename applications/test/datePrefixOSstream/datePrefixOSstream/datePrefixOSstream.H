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

Class
    Foam::datePrefixOSstream

Description
    Version of OSstream which prints a prefix on each line.

    This is useful for running in parallel as it allows the processor number
    to be automatically prepended to each message line.

SourceFiles
    datePrefixOSstream.C

\*---------------------------------------------------------------------------*/

#ifndef datePrefixOSstream_H
#define datePrefixOSstream_H

#include "prefixOSstream2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class datePrefixOSstream Declaration
\*---------------------------------------------------------------------------*/

class datePrefixOSstream
:
    public prefixOSstream2
{
protected:

    // Private Member Functions

        virtual void writePrefix();


public:

    // Constructors

        //- Set stream status
        datePrefixOSstream
        (
            ostream& os,
            const string& name,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Global predefined stream

extern datePrefixOSstream Pdate;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


//- Report a variable name and value
//  using Foam::Pdate in file __FILE__ at line __LINE__
#define DateVar(var)                                                          \
{                                                                             \
    ::Foam::string oldPrefix(::Foam::Pdate.prefix());                         \
    ::Foam::Pdate<< "["<< __FILE__ << ":" << __LINE__ << "] " #var "=";       \
    ::Foam::Pdate.prefix() = oldPrefix + #var " ";                            \
    ::Foam::Pdate<< var << ::Foam::endl;                                      \
    ::Foam::Pdate.prefix() = oldPrefix;                                       \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
