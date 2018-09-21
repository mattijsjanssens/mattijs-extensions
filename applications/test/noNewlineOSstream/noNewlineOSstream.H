/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::noNewlineOSstream

Description
    Version of OSstream which does not print newlines.

SourceFiles
    noNewlineOSstream.C

\*---------------------------------------------------------------------------*/

#ifndef noNewlineOSstream_H
#define noNewlineOSstream_H

#include "OSstream.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class noNewlineOSstream Declaration
\*---------------------------------------------------------------------------*/

class noNewlineOSstream
:
    public OSstream
{
public:

    // Constructors

        //- Set stream status
        noNewlineOSstream
        (
            OSstream& os,
            const string& name = "noNewlineOSstream",
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    // Member functions

        // Write functions

            //- Write character
            virtual Ostream& write(const char);


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream&) const;


    // Member operators

        //- Ostream Operator
        template<class T>
        friend Ostream& operator<<
        (
            const noNewlineOSstream& os,
            const T& l
        )
        {
            return Foam::operator<<(os(), l);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
