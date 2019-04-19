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
    OSstream& os_;


public:

    // Constructors

        //- Set stream status
        explicit noNewlineOSstream
        (
            OSstream& os,
            const string& name = "noNewlineOSstream",
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        )
        :
            OSstream(os.stdStream(), name, format, version, compression),
            os_(os)
        {}


    // Member functions

        // Enquiry

            //- Return the name of the stream
            //  Useful for Fstream to return the filename
            virtual const fileName& name() const
            {
                return os_.name();
            }

            //- Return non-const access to the name of the stream
            //  Useful to alter the stream name
            virtual fileName& name()
            {
                return os_.name();
            }

            //- Return flags of output stream
            virtual ios_base::fmtflags flags() const
            {
                return os_.flags();
            }


        // Write functions

            //- Write next token to stream
            virtual Ostream& write(const token& t)
            {
                return os_.write(t);
            }

            //- Write character
            virtual Ostream& write(const char c)
            {
                if (c == token::NL)
                {
                    return os_.write(' ');
                }
                else
                {
                    return os_.write(c);
                }
            }

            //- Write character string
            virtual Ostream& write(const char* s)
            {
                return os_.write(s);
            }

            //- Write word
            virtual Ostream& write(const word& w)
            {
                return os_.write(w);
            }

            //- Write string
            virtual Ostream& write(const string& s)
            {
                return os_.write(s);
            }

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            virtual Ostream& writeQuoted
            (
                const std::string& s,
                const bool quoted=true
            )
            {
                return os_.writeQuoted(s, quoted);
            }

            //- Write int32_t
            virtual Ostream& write(const int32_t v)
            {
                return os_.write(v);
            }

            //- Write int64_t
            virtual Ostream& write(const int64_t v)
            {
                return os_.write(v);
            }

            //- Write floatScalar
            virtual Ostream& write(const floatScalar v)
            {
                return os_.write(v);
            }

            //- Write doubleScalar
            virtual Ostream& write(const doubleScalar v)
            {
                return os_.write(v);
            }

            //- Write longDoubleScalar
            virtual Ostream& write(const longDoubleScalar v)
            {
                return os_.write(v);
            }

            //- Write binary block
            virtual Ostream& write(const char* s, std::streamsize sz)
            {
                return os_.write(s, sz);
            }

            //- Add indentation characters
            virtual void indent()
            {
                os_.indent();
            }


        // Stream state functions

            //- Set flags of output stream
            virtual ios_base::fmtflags flags(const ios_base::fmtflags flags)
            {
                return os_.flags(flags);
            }

            //- Flush stream
            virtual void flush()
            {
                os_.flush();
            }

            //- Add newline and flush stream
            virtual void endl()
            {
                os_.endl();
            }

            //- Get width of output field
            virtual int width() const
            {
                return os_.width();
            }

            //- Set width of output field (and return old width)
            virtual int width(const int n)
            {
                return os_.width(n);
            }

            //- Get precision of output field
            virtual int precision() const
            {
                return os_.precision();
            }

            //- Set precision of output field (and return old precision)
            virtual int precision(const int n)
            {
                return os_.precision(n);
            }


        // STL stream

            //- Access to underlying std::ostream
            virtual ostream& stdStream()
            {
                return os_.stdStream();
            }

            //- Const access to underlying std::ostream
            virtual const ostream& stdStream() const
            {
                return os_.stdStream();
            }


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream& os) const
            {
                os  << "noNewlineOSstream ";
                OSstream::print(os);
            }


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
