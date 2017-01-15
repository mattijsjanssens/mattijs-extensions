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
    Test-volField

Description
    The idea is that argList allocates a fileServer (through command-line
    arguments or environment vars). This has operations to do
    - file existence checking
    - mkDir, rm etc.
    - open an IFstream
    - open an OFstream
    and everywhere instead of constructing an I/OFstream we do a
        fileServer::NewIFstream(io)
        fileServer::NewOFstream(io)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IOPtrList.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "decomposedBlockData.H"

using namespace Foam;

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOList<List<char>>, 0);
}

template<class T>
bool read(autoPtr<Istream>& isPtr, T& localData)
{
    bool ok = false;

    if (Pstream::master())
    {
        Istream& is = isPtr();

        is.fatalCheck("read(Istream&)");

        token firstToken(is);

        is.fatalCheck
        (
            "read(Istream&) : "
            "reading first token"
        );

        if (firstToken.isLabel())
        {
            // Read size of list
            label s = firstToken.labelToken();

            if (s != Pstream::nProcs())
            {
                FatalIOErrorInFunction
                (
                    is
                )   << "size " << s
                    << " differs from the number of processors "
                    << Pstream::nProcs()
                    << firstToken.info()
                    << exit(FatalIOError);
            }

            // Read beginning of contents
            char delimiter = is.readBeginList("List");

            if (s)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    // Read master data
                    {
                        Pout<< "Reading data for master" << endl;
                        is >> localData;
                        is.fatalCheck
                        (
                            "read(Istream&) : "
                            "reading entry"
                        );
                        Pout<< "Finished reading data for master" << endl;
                    }

                    for (label proci = 1; proci < s; proci++)
                    {
                        Pout<< "Reading file for processor " << proci << endl;

                        T elems(is);

                        Pout<< "Read file data:" << elems << endl;

                        is.fatalCheck
                        (
                            "read(Istream&) : "
                            "reading entry"
                        );

                        OPstream os(Pstream::scheduled, proci);
                        os << elems;
                    }
                }
                else
                {
                    FatalIOErrorInFunction
                    (
                        is
                    )   << "incorrect delimiter, expected (, found "
                        << delimiter
                        << exit(FatalIOError);
                }
            }

            // Read end of contents
            is.readEndList("List");

            ok = is.good();
        }
        else
        {
            FatalIOErrorInFunction
            (
                is
            )   << "incorrect first token, expected <int>, found "
                << firstToken.info()
                << exit(FatalIOError);
        }
    }
    else
    {
        Pout<< "Receiving data from master" << endl;

        IPstream is(Pstream::scheduled, Pstream::masterNo());
        is >> localData;
        Pout<< "Received data:" << localData << endl;
    }

    Pstream::scatter(ok);

    return ok;
}


bool write
(
    autoPtr<OSstream>& osPtr,
    const UList<char>& L,
    List<std::streamoff>& start
)
{
    bool ok = false;

    if (Pstream::master())
    {
        start.setSize(Pstream::nProcs());

        OSstream& os = osPtr();

        // Write size and start delimiter
        os << nl << Pstream::nProcs() << nl << token::BEGIN_LIST;

        // Write master data
        {
            // os << nl << L.size() << nl;
            // if (L.size())
            // {
            //     os.write
            //     (
            //         reinterpret_cast<const char*>(L.begin()),
            //         L.byteSize()
            //     );
            // }
            os << nl;
            start[Pstream::masterNo()] = os.stdStream().tellp();
            os << L;
        }
        // Write slaves
        for (label proci = 1; proci < Pstream::nProcs(); proci++)
        {
            Pout<< "Receiving and writing data for processor " << proci << endl;

            IPstream is(Pstream::scheduled, proci);
            List<char> elems(is);

            is.fatalCheck
            (
                "write(Istream&) : "
                "reading entry"
            );

            // os << nl << elems.size() << nl;
            // if (elems.size())
            // {
            //     os.write
            //     (
            //         reinterpret_cast<const char*>(elems.begin()),
            //         elems.byteSize()
            //     );
            // }
            os << nl;
            start[proci] = os.stdStream().tellp();
            os << elems;
        }

        // Write end delimiter
        os << nl << token::END_LIST << nl;

        ok = os.good();
    }
    else
    {
        Pout<< "Sending data to processor " << Pstream::masterNo() << endl;

        OPstream os(Pstream::scheduled, Pstream::masterNo());
        os << L;

        Pout<< "Done Sending data to processor " << Pstream::masterNo() << endl;
    }

Pout<< "scatter:" << ok << endl;
    Pstream::scatter(ok);

    return ok;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

Pout<< "std::streamoff:" << sizeof(std::streamoff) << endl;
Pout<< "off_t:" << sizeof(off_t) << endl;
Pout<< "label:" << sizeof(label) << endl;


//     IOobject io
//     (
//         "bufs",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE,
//         false
//     );
//
//     {
//         // Allocate
//         IOList<List<char>> bufs(io, 2);
//
//         // Fill
//         forAll(bufs, proci)
//         {
//             //bufs.set(proci, new List<char>(0));
//
//             List<char>& buf = bufs[proci];
//             buf.setSize(10);
//             forAll(buf, i)
//             {
//                 buf[i] = proci+'a'+i;
//             }
//         }
//
//         // Write
//         bufs.write();
//     }



    IOobject io
    (
        "bufs2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    // Write a field
    {
        volScalarField fld
        (
            IOobject
            (
                "myFld",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );
        fld == dimensionedScalar("zero", dimless, 0.0);

        // Test normal writing
        {
            fld.write();
        }


        // Test collating writing
        OStringStream fldStream;
        fldStream << fld;

        string s(fldStream.str());
        //UList<char> localData(const_cast<char*>(s.data()), s.size());
        //autoPtr<OSstream> osPtr;
        //if (Pstream::master())
        //{
        //    osPtr.reset(new OFstream(io.objectPath(), IOstream::BINARY));
        //    io.writeHeader(osPtr());
        //}
        //List<std::streamoff> start;
        //write(osPtr, localData, start);
        //DebugVar(start);
        UList<char> slice(const_cast<char*>(s.data()), label(s.size()));
        List<char> sList(slice);
        decomposedBlockData localData(io, sList.xfer());
        localData.write();
    }
    {
        io.readOpt() = IOobject::MUST_READ;
        decomposedBlockData localData(io);

        string s(localData.begin(), localData.size());
        IStringStream is(s);
        dictionary dict(is);
        volScalarField fld
        (
            IOobject
            (
                "myFld2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dict
        );
        DebugVar(fld);
    }


//
//
//     // Scan for start and size of entries
//     List<char> localData;
//     {
//         IOList<List<char>> bufs(io);
//
//         autoPtr<Istream> isPtr;
//         if (Pstream::master())
//         {
//             isPtr.reset(new IFstream(io.objectPath()));
//             io.readHeader(isPtr());
//         }
//         read(isPtr, localData);
//
//         DebugVar(localData);
//     }
//
//     // Read a dictionary from the stream
//     {
//         string s(localData.begin(), localData.size());
// DebugVar(s);
//         IStringStream is(s);
//         dictionary dict(is);
//         DebugVar(dict);
//         volScalarField fld
//         (
//             IOobject
//             (
//                 "myFld2",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE,
//                 false
//             ),
//             mesh,
//             dict
//         );
//     }

    Pout<< "**end**" << endl;


    return 0;
}


// ************************************************************************* //
