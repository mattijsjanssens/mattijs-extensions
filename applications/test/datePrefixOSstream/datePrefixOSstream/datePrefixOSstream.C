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

#include "datePrefixOSstream.H"
#include <sys/time.h>

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::datePrefixOSstream::writePrefix()
{
    // Loosely based on
    // http://stackoverflow.com/questions/8304259/formatting-struct-timespec

    const uint TIME_FMT = strlen("2012-12-31 12:59:59.123456789") + 1;
    char timestr[TIME_FMT];

    // Get seconds since 1970
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);

    // Convert into day, hours, second etc.
    struct tm t;
    localtime_r(&(ts.tv_sec), &t);

    // Format as string
    int nChar = strftime(timestr, TIME_FMT, "%F %T", &t);
    // Append nanoseconds
    snprintf(&timestr[nChar], TIME_FMT-nChar, ".%09ld", ts.tv_nsec);

    if (prefix().size())
    {
        OSstream::write(prefix().c_str());
    }
    OSstream::write('[');
    OSstream::write(timestr);
    OSstream::write("] ");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::datePrefixOSstream::datePrefixOSstream
(
    ostream& os,
    const string& name,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    prefixOSstream2(os, name, format, version, compression)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::datePrefixOSstream::print(Ostream& os) const
{
    os  << "datePrefixOSstream ";
    prefixOSstream2::print(os);
}


/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    datePrefixOSstream Pdate(cout, "Pdate");
}


// ************************************************************************* //
