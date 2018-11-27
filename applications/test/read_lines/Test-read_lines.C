/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Mattijs Janssens
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

Description:
    See
    https://stackoverflow.com/questions/7868936/
        read-file-line-by-line-using-ifstream-in-c

    See also performance table in
https://www.boost.org/doc/libs/1_47_0/libs/conversion/lexical_cast.htm#tuning

    Written OpenFOAM 40000000 in 15
    Written C (vectorised) 40000000 in 8.4
    Written C 40000000 in 9.42
    Read C++ (line based) 40000000 in 12.22
    Read OpenFOAM (stream) 40000000 in 10.12
    Read C++ (stream) 40000000 in 13.18
    Read C (line based) 40000000 in 10.71



\*---------------------------------------------------------------------------*/

#include "UList.H"
#include "List.H"
#include "fvMesh.H"
#include "volFields.H"
#include "argList.H"
#include "Time.H"
#include "OFstream.H"
#include "IFstream.H"
#include "Random.H"
#include "IOmanip.H"
#include <boost/lexical_cast.hpp>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
void Foam::UList<double>::writeNonUniform(Ostream& os) const
{
//Pout<< "width:" << os.width() << endl;
//Pout<< "precision:" << os.precision() << endl;
//Pout<< "fmtflags:" << os.flags() << endl;
//Pout<< "boolalpha:" << (os.flags() & ios_base::boolalpha) << endl;
//Pout<< "showbase:" << (os.flags() & ios_base::showbase) << endl;
//Pout<< "showpoint:" << (os.flags() & ios_base::showpoint) << endl;
//Pout<< "showpos:" << (os.flags() & ios_base::showpos) << endl;
//Pout<< "skipws:" << (os.flags() & ios_base::skipws) << endl;
//Pout<< "unitbuf:" << (os.flags() & ios_base::unitbuf) << endl;
//Pout<< "uppercase:" << (os.flags() & ios_base::uppercase) << endl;
//
//Pout<< "dec:" << (os.flags() & ios_base::showbase) << endl;
//Pout<< "hex:" << (os.flags() & ios_base::hex) << endl;
//Pout<< "oct:" << (os.flags() & ios_base::oct) << endl;
//
//Pout<< "fixed:" << (os.flags() & ios_base::fixed) << endl;
//Pout<< "scientific:" << (os.flags() & ios_base::scientific) << endl;
//
//Pout<< "internal:" << (os.flags() & ios_base::internal) << endl;
//Pout<< "left:" << (os.flags() & ios_base::left) << endl;
//Pout<< "right:" << (os.flags() & ios_base::right) << endl;

    const Foam::UList<double>& L = *this;

    // Write size and start delimiter
    os << nl << L.size() << nl << token::BEGIN_LIST;

    if
    (
        (os.width() == 0)
     && !(os.flags() & ios_base::showbase)
     && !(os.flags() & ios_base::showpoint)
     && !(os.flags() & ios_base::showpos)
     && !(os.flags() & ios_base::unitbuf)
     && !(os.flags() & ios_base::uppercase)
     && !(os.flags() & ios_base::hex)
     && !(os.flags() & ios_base::oct)
     && !(os.flags() & ios_base::fixed)
     && !(os.flags() & ios_base::scientific)
     && !(os.flags() & ios_base::internal)
     && !(os.flags() & ios_base::right)
     && isA<OFstream>(os)
    )
    {
        OFstream& ofs = dynamic_cast<OFstream&>(os);

        // Can do optimised
        Pout<< "Optimised" << endl;
        const string singleFmt("\n%-." + Foam::name(ofs.precision()) + "lg");
        string fmt;
        for (int i = 0; i < 8; i++)
        {
            fmt += singleFmt;
        }
        //Pout<< "Format string:" << fmt << endl;
        const char* fmt_ptr = fmt.c_str();

        const int bufLen = 4096;
        char buf[bufLen];
        int free = 0;
        label i = 0;
        for (label outer = 0; outer < L.size()/8; outer++)
        {
            const int nFree = bufLen-free;
            int n = snprintf
            (
                &buf[free],
                nFree,
                fmt_ptr,
                L[i],
                L[i+1],
                L[i+2],
                L[i+3],
                L[i+4],
                L[i+5],
                L[i+6],
                L[i+7]
            );
            if (n >= nFree)
            {
                // Failed. Flush buffer so far
                ofs.stdStream().write(buf, free);
                free = 0;

                // Re-print current entry (assume it succeeds)
                n = snprintf
                (
                    &buf[free],
                    bufLen-free,
                    fmt_ptr,
                    L[i],
                    L[i+1],
                    L[i+2],
                    L[i+3],
                    L[i+4],
                    L[i+5],
                    L[i+6],
                    L[i+7]
                );
            }
            i += 8;
            free += n;
        }

        // Do remainder element by element
        const char* sptr = singleFmt.c_str();
        while (i < L.size())
        {
            const int nFree = bufLen-free;
            int n = snprintf(&buf[free], nFree, sptr, L[i]);
            if (n >= nFree)
            {
                ofs.stdStream().write(buf, free);
                free = 0;

                // Re-print current entry
                n = snprintf(&buf[free], bufLen-free, sptr, L[i]);
            }
            free += n;
            i++;
        }

        if (free > 0)
        {
            ofs.stdStream().write(buf, free);
        }
    }
    else
    {
        // Write contents
        forAll(L, i)
        {
            os << nl << L[i];
        }
    }

    // Write end delimiter
    os << nl << token::END_LIST << nl;
}


template<>
void Foam::UList<Foam::Vector<double>>::writeNonUniform(Ostream& os) const
{
//Pout<< "width:" << os.width() << endl;
//Pout<< "precision:" << os.precision() << endl;
//Pout<< "fmtflags:" << os.flags() << endl;
//Pout<< "boolalpha:" << (os.flags() & ios_base::boolalpha) << endl;
//Pout<< "showbase:" << (os.flags() & ios_base::showbase) << endl;
//Pout<< "showpoint:" << (os.flags() & ios_base::showpoint) << endl;
//Pout<< "showpos:" << (os.flags() & ios_base::showpos) << endl;
//Pout<< "skipws:" << (os.flags() & ios_base::skipws) << endl;
//Pout<< "unitbuf:" << (os.flags() & ios_base::unitbuf) << endl;
//Pout<< "uppercase:" << (os.flags() & ios_base::uppercase) << endl;
//
//Pout<< "dec:" << (os.flags() & ios_base::showbase) << endl;
//Pout<< "hex:" << (os.flags() & ios_base::hex) << endl;
//Pout<< "oct:" << (os.flags() & ios_base::oct) << endl;
//
//Pout<< "fixed:" << (os.flags() & ios_base::fixed) << endl;
//Pout<< "scientific:" << (os.flags() & ios_base::scientific) << endl;
//
//Pout<< "internal:" << (os.flags() & ios_base::internal) << endl;
//Pout<< "left:" << (os.flags() & ios_base::left) << endl;
//Pout<< "right:" << (os.flags() & ios_base::right) << endl;

    const Foam::UList<Vector<double>>& L = *this;

    // Write size and start delimiter
    os << nl << L.size() << nl << token::BEGIN_LIST;

    if
    (
        (os.width() == 0)
     && !(os.flags() & ios_base::showbase)
     && !(os.flags() & ios_base::showpoint)
     && !(os.flags() & ios_base::showpos)
     && !(os.flags() & ios_base::unitbuf)
     && !(os.flags() & ios_base::uppercase)
     && !(os.flags() & ios_base::hex)
     && !(os.flags() & ios_base::oct)
     && !(os.flags() & ios_base::fixed)
     && !(os.flags() & ios_base::scientific)
     && !(os.flags() & ios_base::internal)
     && !(os.flags() & ios_base::right)
     && isA<OFstream>(os)
    )
    {
        OFstream& ofs = dynamic_cast<OFstream&>(os);

        // Can do optimised
        Pout<< "Optimised" << endl;
        const string sFmt("%-." + Foam::name(ofs.precision()) + "lg");
        string fmt("\n("+sFmt+' '+sFmt+' '+sFmt+')');
        Pout<< "Format string:" << fmt << endl;
        const char* fmt_ptr = fmt.c_str();

        const int bufLen = 4096;
        char buf[bufLen];
        int free = 0;
        forAll(L, i)
        {
            const Vector<double>& p = L[i];

            const int nFree = bufLen-free;
            int n = snprintf(&buf[free], nFree, fmt_ptr, p.x(), p.y(), p.z());
            if (n >= nFree)
            {
                ofs.stdStream().write(buf, free);
                free = 0;

                // Re-print current entry
                n = snprintf(&buf[free], bufLen, fmt_ptr, p.x(), p.y(), p.z());
            }
            free += n;
        }

        if (free > 0)
        {
            ofs.stdStream().write(buf, free);
        }
    }
    else
    {
        // Write contents
        forAll(L, i)
        {
            os << nl << L[i];
        }
    }

    // Write end delimiter
    os << nl << token::END_LIST << nl;
}


// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    //#include "createMesh.H"

    Random rndGen(0);

    scalarField fld(10000000);
    forAll(fld, i)
    {
        fld[i] = pow(10, 10*rndGen.scalar01());
    }
    vectorField vfld(4000000);
    forAll(vfld, i)
    {
        vfld[i] = rndGen.sample01<vector>();
    }

    // Write scalarField in ascii
    cpuTime timer;

    {
        OFstream os("scalars_foam.txt");
        os << fld << endl;
    }
    Pout<< "Written OpenFOAM " << fld.size() << " in "
        << timer.cpuTimeIncrement()
        << endl;

    {
        OFstream os("scalars_snprintf.txt");
        os.precision(IOstream::defaultPrecision());

        const string fmt("%-." + Foam::name(os.precision()) + "lg\n");
        //Pout<< "Format string:" << fmt << endl;
        const char* fmt_ptr = fmt.c_str();

        const int bufLen = 4096;
        char buf[bufLen];
        int free = 0;
        forAll(fld, i)
        {
            const int nFree = bufLen-free;
            int n = snprintf(&buf[free], nFree, fmt_ptr, fld[i]);
            if (n >= nFree)
            {
                os.stdStream().write(buf, free);
                free = 0;

                // Re-print current entry
                n = snprintf(&buf[free], 20, fmt_ptr, fld[i]);
            }
            free += n;
        }

        if (free > 0)
        {
            os.stdStream().write(buf, free);
        }
    }
    Pout<< "Written snprintf" << fld.size() << " in "
        << timer.cpuTimeIncrement() << endl;

    {
        OFstream os("scalars_snprintf_vectorised.txt");
        os.precision(IOstream::defaultPrecision());

        const string singleFmt("%-." + Foam::name(os.precision()) + "lg\n");
        string fmt;
        for (int i = 0; i < 8; i++)
        {
            fmt += singleFmt;
        }
        //Pout<< "Format string:" << fmt << endl;
        const char* fmt_ptr = fmt.c_str();

        const int bufLen = 4096;
        char buf[bufLen];
        int free = 0;
        label i = 0;
        for (label outer = 0; outer < fld.size()/8; outer++)
        {
            const int nFree = bufLen-free;
            int n = snprintf
            (
                &buf[free],
                nFree,
                fmt_ptr,
                fld[i],
                fld[i+1],
                fld[i+2],
                fld[i+3],
                fld[i+4],
                fld[i+5],
                fld[i+6],
                fld[i+7]
            );
            if (n >= nFree)
            {
                // Failed. Flush buffer so far
                os.stdStream().write(buf, free);
                free = 0;

                // Re-print current entry (assume it succeeds)
                n = snprintf
                (
                    &buf[free],
                    bufLen-free,
                    fmt_ptr,
                    fld[i],
                    fld[i+1],
                    fld[i+2],
                    fld[i+3],
                    fld[i+4],
                    fld[i+5],
                    fld[i+6],
                    fld[i+7]
                );
            }
            i += 8;
            free += n;
        }

        const char* sptr = singleFmt.c_str();
        while (i < fld.size())
        {
            const int nFree = bufLen-free;
            int n = snprintf(&buf[free], nFree, sptr, fld[i]);
            if (n >= nFree)
            {
                os.stdStream().write(buf, free);
                free = 0;

                // Re-print current entry
                n = snprintf(&buf[free], bufLen-free, sptr, fld[i]);
            }
            free += n;
            i++;
        }

        if (free > 0)
        {
            os.stdStream().write(buf, free);
        }
    }
    Pout<< "Written snprintf (vectorised)" << fld.size() << " in "
        << timer.cpuTimeIncrement() << endl;




//- OStringStream is slower (7s vs. 4s)
//    {
//        OStringStream buf;
//        buf << fld << endl;
//        OFstream os("scalars_foam.txt");
//        os << buf.str() << endl;
//    }
//    Pout<< "Written buffered OpenFOAM " << fld.size() << " in "
//        << timer.cpuTimeIncrement()
//        << endl;

    {
        FILE* fp = fopen("scalars_c.txt", "w");
        const string singleFmt
        (
            "\n%-." + Foam::name(IOstream::defaultPrecision()) + "lg"
        );
        const char* fmt_ptr = singleFmt.c_str();
        forAll(fld, i)
        {
            fprintf(fp, fmt_ptr, fld[i]);
        }
        fclose(fp);
    }
    Pout<< "Written C " << fld.size() << " in " << timer.cpuTimeIncrement()
        << endl;

    {
        FILE* fp = fopen("scalars_c.txt", "w");

        const string singleFmt
        (
            "\n%-." + Foam::name(IOstream::defaultPrecision()) + "lg"
        );
        string fmt;
        for (int i = 0; i < 8; i++)
        {
            fmt += singleFmt;
        }
        //Pout<< "Format string:" << fmt << endl;
        const char* fmt_ptr = fmt.c_str();

        label i = 0;
        for (label outer = 0; outer < fld.size()/8; outer++)
        {
            //Pout<< "Outer:" << outer << " i:" << i << endl;
            fprintf
            (
                fp,
                fmt_ptr,
                fld[i],
                fld[i+1],
                fld[i+2],
                fld[i+3],
                fld[i+4],
                fld[i+5],
                fld[i+6],
                fld[i+7]
            );

            i += 8;
        }
        const char* sptr = singleFmt.c_str();
        while (i < fld.size())
        {
            fprintf(fp, sptr, fld[i++]);
        }
        fclose(fp);
    }
    Pout<< "Written C (vectorised) " << fld.size() << " in "
        << timer.cpuTimeIncrement() << endl;

    // Read
    {
        std::ifstream file("scalars_c.txt");
        if (file.is_open())
        {
            label i = 0;
            std::string line;
            while (getline(file, line))
            {
                sscanf(line.c_str(), "%lg", &fld[i]);
                i++;
            }
            file.close();
        }
    } 
    Pout<< "Read C++ (line based) " << fld.size() << " in "
        << timer.cpuTimeIncrement() << endl;

    {
        IFstream file("scalars_c.txt");
        if (file.good())
        {
            label i = 0;
            while (i < fld.size())
            {
                file >> fld[i];
                i++;
            }
        }
    } 
    Pout<< "Read OpenFOAM (stream) " << fld.size() << " in "
        << timer.cpuTimeIncrement() << endl;

    {
        std::ifstream file("scalars_c.txt");
        if (file.good())
        {
            label i = 0;
            while (i < fld.size())
            {
                file >> fld[i];
                i++;
            }
        }
    } 
    Pout<< "Read C++ (stream) " << fld.size()
        << " in " << timer.cpuTimeIncrement() << endl;

//    {
//        std::ifstream file("scalars_c.txt");
//        if (file.good())
//        {
//            // get length of file:
//            file.seekg (0, file.end);
//            std::streamoff length = file.tellg();
//            file.seekg (0, file.beg);
//
//            char * buffer = new char [length];
//            // read data as a block:
//            file.read (buffer,length);
//
//            Pout<< "Read file (C++) " << length << " in "
//                << timer.cpuTimeIncrement() << endl;
//
//            IStringStream is(buffer);
//
//            label i = 0;
//            while (is.good())
//            {
//                is >> fld[i];
//                i++;
//            }
//        }
//        Pout<< "Parsed OpenFOAM " << fld.size() << " in "
//            << timer.cpuTimeIncrement()
//            << endl;
//    }

    {
        FILE* fp = fopen("scalars_c.txt", "r");
        if (fp == NULL)
            exit(1);

        char* line = NULL;
        size_t len = 0;

        label i = 0;
        while ((getline(&line, &len, fp)) != -1) {
            // using printf() in all tests for consistency

            sscanf(line, "%lg", &fld[i]);
            i++;
        }
        fclose(fp);
        if (line)
            free(line);
    }
    Pout<< "Read C (line based) " << fld.size() << " in "
        << timer.cpuTimeIncrement() << endl;


    // Vectors

    {
        OFstream os("vectors_foam.txt");
        os << vfld << endl;
    }
    Pout<< "VECTOR Written OpenFOAM " << fld.size() << " in "
        << timer.cpuTimeIncrement()
        << endl;

    {
        FILE* fp = fopen("vectors_c.txt", "w");

        const string sFmt
        (
            "%-." + Foam::name(IOstream::defaultPrecision()) + "lg"
        );
        string fmt("\n("+sFmt+' '+sFmt+' '+sFmt+')');
        Pout<< "Format string:" << fmt << endl;
        const char* fmt_ptr = fmt.c_str();

        forAll(vfld, i)
        {
            const vector& p = vfld[i];
            fprintf(fp, fmt_ptr, p.x(), p.y(), p.z());
        }
        fclose(fp);
    }
    Pout<< "VECTOR Written C " << fld.size() << " in "
        << timer.cpuTimeIncrement() << endl;



    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
