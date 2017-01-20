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

#include "decomposedBlockData.H"
#include "Pstream.H"
#include "IOstreams.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamBuffers.H"
#include "OFstream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decomposedBlockData, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decomposedBlockData::decomposedBlockData
(
    const IOobject& io,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "IOList " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOList does not support automatic rereading."
            << endl;
    }
    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read();
    }
}


Foam::decomposedBlockData::decomposedBlockData
(
    const IOobject& io,
    const List<char>& list,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read();
    }
    else
    {
        List<char>::operator=(list);
    }
}


Foam::decomposedBlockData::decomposedBlockData
(
    const IOobject& io,
    const Xfer<List<char>>& list,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "IOList " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOList does not support automatic rereading."
            << endl;
    }

    List<char>::transfer(list());

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read();
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::decomposedBlockData::~decomposedBlockData()
{}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

bool Foam::decomposedBlockData::readBlocks
(
    autoPtr<Istream>& isPtr,
    List<char>& data,
    const UPstream::commsTypes commsType
)
{
    bool ok = false;

    if (commsType == UPstream::scheduled)
    {
        if (UPstream::master())
        {
            Istream& is = isPtr();
            is.fatalCheck("read(Istream&)");
            token firstToken(is);
            is.fatalCheck("read(Istream&) : " "reading first token");

            if (firstToken.isLabel())
            {
                // Read size of list
                label s = firstToken.labelToken();

                if (s != UPstream::nProcs())
                {
                    FatalIOErrorInFunction
                    (
                        is
                    )   << "size " << s
                        << " differs from the number of processors "
                        << UPstream::nProcs() << exit(FatalIOError);
                }

                // Read beginning of contents
                char delimiter = is.readBeginList("List");

                if (s)
                {
                    if (delimiter == token::BEGIN_LIST)
                    {
                        // Read master data
                        {
                            is >> data;
                            is.fatalCheck
                            (
                                "read(Istream&) : "
                                "reading entry"
                            );
                        }

                        for (label proci = 1; proci < s; proci++)
                        {
                            List<char> elems(is);
                            is.fatalCheck
                            (
                                "read(Istream&) : "
                                "reading entry"
                            );

                            OPstream os(UPstream::scheduled, proci);
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
                    << firstToken.info() << exit(FatalIOError);
            }
        }
        else
        {
            IPstream is(UPstream::scheduled, UPstream::masterNo());
            is >> data;
        }
    }
    else
    {
        PstreamBuffers pBufs(UPstream::nonBlocking);

        if (UPstream::master())
        {
            Istream& is = isPtr();
            is.fatalCheck("read(Istream&)");
            token firstToken(is);
            is.fatalCheck("read(Istream&) : " "reading first token");

            if (firstToken.isLabel())
            {
                // Read size of list
                label s = firstToken.labelToken();

                if (s != UPstream::nProcs())
                {
                    FatalIOErrorInFunction
                    (
                        is
                    )   << "size " << s
                        << " differs from the number of processors "
                        << UPstream::nProcs() << exit(FatalIOError);
                }

                // Read beginning of contents
                char delimiter = is.readBeginList("List");

                if (s)
                {
                    if (delimiter == token::BEGIN_LIST)
                    {
                        // Read master data
                        {
                            is >> data;
                            is.fatalCheck
                            (
                                "read(Istream&) : "
                                "reading entry"
                            );
                        }

                        for (label proci = 1; proci < s; proci++)
                        {
                            List<char> elems(is);

                            is.fatalCheck
                            (
                                "read(Istream&) : "
                                "reading entry"
                            );

                            UOPstream os(proci, pBufs);
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
                    << firstToken.info() << exit(FatalIOError);
            }
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (!UPstream::master())
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> data;
        }
    }

    Pstream::scatter(ok);

    return ok;
}


bool Foam::decomposedBlockData::writeBlocks
(
    autoPtr<OSstream>& osPtr,
    List<std::streamoff>& start,
    const UList<char>& data,
    const UPstream::commsTypes commsType
)
{
    bool ok = false;

    if (commsType == UPstream::scheduled)
    {
        if (UPstream::master())
        {
            start.setSize(UPstream::nProcs());

            OSstream& os = osPtr();

            // Write size and start delimiter
            os << nl << UPstream::nProcs() << nl << token::BEGIN_LIST;

            // Write master data
            {
                os << nl;
                start[UPstream::masterNo()] = os.stdStream().tellp();
                os << data;
            }
            // Write slaves
            for (label proci = 1; proci < UPstream::nProcs(); proci++)
            {
                IPstream is(UPstream::scheduled, proci);
                List<char> elems(is);

                is.fatalCheck("write(Istream&) : reading entry");

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
            OPstream os(UPstream::scheduled, UPstream::masterNo());
            os << data;
        }
    }
    else
    {
        PstreamBuffers pBufs(UPstream::nonBlocking);

        if (!UPstream::master())
        {
            UOPstream os(UPstream::masterNo(), pBufs);
            os << data;
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (UPstream::master())
        {
            start.setSize(UPstream::nProcs());

            OSstream& os = osPtr();

            // Write size and start delimiter
            os << nl << UPstream::nProcs() << nl << token::BEGIN_LIST;

            // Write master data
            {
                os << nl;
                start[UPstream::masterNo()] = os.stdStream().tellp();
                os << data;
            }
            // Write slaves
            for (label proci = 1; proci < UPstream::nProcs(); proci++)
            {
                UIPstream is(UPstream::masterNo(), pBufs);
                List<char> elems(is);

                is.fatalCheck("write(Istream&) : reading entry");

                os << nl;
                start[proci] = os.stdStream().tellp();
                os << elems;
            }

            // Write end delimiter
            os << nl << token::END_LIST << nl;

            ok = os.good();
        }
    }

    Pstream::scatter(ok);

    return ok;
}


bool Foam::decomposedBlockData::read()
{
    autoPtr<Istream> isPtr;
    if (UPstream::master())
    {
        isPtr.reset(new IFstream(objectPath()));
        readHeader(isPtr());
    }
    return readBlocks(isPtr, *this, commsType_);
}


bool Foam::decomposedBlockData::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
      autoPtr<OSstream> osPtr;
      if (UPstream::master())
      {
          // Note: always write binary. These are strings so readable
          //       anyway. They have already be tokenised on the sending side.
          osPtr.reset(new OFstream(objectPath(), IOstream::BINARY));
          writeHeader(osPtr());
      }
      List<std::streamoff> start;
      return writeBlocks(osPtr, start, *this, commsType_);
}


// ************************************************************************* //
