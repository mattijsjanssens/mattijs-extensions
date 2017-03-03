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
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamBuffers.H"
#include "OFstream.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decomposedBlockData, 0);
}

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
}


Foam::decomposedBlockData::decomposedBlockData
(
    const IOobject& io,
    const UList<char>& list,
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
            << "decomposedBlockData " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
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

bool Foam::decomposedBlockData::readMasterHeader(IOobject& io, Istream& is)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readMasterHeader:"
            << " stream:" << is.name() << endl;
    }

    // Master-only reading of header
    is.fatalCheck("read(Istream&)");

    List<char> data(is);
    is.fatalCheck("read(Istream&) : reading entry");
    string buf(data.begin(), data.size());
    IStringStream str(is.name(), buf);

    return io.readHeader(str);
}


bool Foam::decomposedBlockData::readBlock
(
    const label blocki,
    Istream& is,
    List<char>& data
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readMasterHeader:"
            << " stream:" << is.name() << " attempt to read block " << blocki
            << endl;
    }

    is.fatalCheck("read(Istream&)");

    for (label i = 0; i < blocki+1; i++)
    {
        // Read data, override old data
        is >> data;
        is.fatalCheck("read(Istream&) : reading entry");
    }
    return is.good();
}


bool Foam::decomposedBlockData::readBlocks
(
    autoPtr<ISstream>& isPtr,
    List<char>& data,
    const UPstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr.valid() ? isPtr().name() : "invalid")
            << " commsType:" << Pstream::commsTypeNames[commsType] << endl;
    }

    bool ok = false;

    if (commsType == UPstream::scheduled)
    {
        if (UPstream::master())
        {
            Istream& is = isPtr();
            is.fatalCheck("read(Istream&)");

            // Read master data
            {
                is >> data;
                is.fatalCheck("read(Istream&) : reading entry");
            }

            for
            (
                label proci = 1;
                proci < UPstream::nProcs();
                proci++
            )
            {
                List<char> elems(is);
                is.fatalCheck("read(Istream&) : reading entry");

                OPstream os(UPstream::scheduled, proci);
                os << elems;
            }

            ok = is.good();
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

            // Read master data
            {
                is >> data;
                is.fatalCheck("read(Istream&) : reading entry");
            }

            for
            (
                label proci = 1;
                proci < UPstream::nProcs();
                proci++
            )
            {
                List<char> elems(is);
                is.fatalCheck("read(Istream&) : reading entry");

                UOPstream os(proci, pBufs);
                os << elems;
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
    if (debug)
    {
        Pout<< "decomposedBlockData::writeBlocks:"
            << " stream:" << (osPtr.valid() ? osPtr().name() : "invalid")
            << " data:" << data.size()
            << " commsType:" << Pstream::commsTypeNames[commsType] << endl;
    }

    bool ok = false;

    if (commsType == UPstream::scheduled)
    {
        if (UPstream::master())
        {
            start.setSize(UPstream::nProcs());

            OSstream& os = osPtr();

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

            // Write master data
            {
                os << nl;
                start[UPstream::masterNo()] = os.stdStream().tellp();
                os << data;
            }

            // Write slaves
            for (label proci = 1; proci < UPstream::nProcs(); proci++)
            {
                UIPstream is(proci, pBufs);
                List<char> elems(is);

                is.fatalCheck("write(Istream&) : reading entry");

                os << nl;
                start[proci] = os.stdStream().tellp();
                os << elems;
            }

            ok = os.good();
        }
    }

    Pstream::scatter(ok);

    return ok;
}


bool Foam::decomposedBlockData::read()
{
    autoPtr<ISstream> isPtr;
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
    IOstream::compressionType cmp,
    const bool valid
) const
{
    autoPtr<OSstream> osPtr;
    if (UPstream::master())
    {
        // Note: always write binary. These are strings so readable
        //       anyway. They have already be tokenised on the sending side.
        osPtr.reset(new OFstream(objectPath(), IOstream::BINARY, ver, cmp));
        writeHeader(osPtr());
    }
    List<std::streamoff> start;
    return writeBlocks(osPtr, start, *this, commsType_);
}


Foam::label Foam::decomposedBlockData::numBlocks(const fileName& fName)
{
    label nBlocks = 0;

    IFstream is(fName);
    is.fatalCheck("decomposedBlockData::numBlocks(const fileName&)");

    if (!is.good())
    {
        return nBlocks;
    }

    // Skip header
    token firstToken(is);

    if
    (
        is.good()
     && firstToken.isWord()
     && firstToken.wordToken() == "FoamFile"
    )
    {
        dictionary headerDict(is);
        is.version(headerDict.lookup("version"));
        is.format(headerDict.lookup("format"));
    }

    List<char> data;
    while (is.good())
    {
        token sizeToken(is);
        if (!sizeToken.isLabel())
        {
            return nBlocks;
        }
        is.putBack(sizeToken);

        is >> data;
        nBlocks++;
    }

    return nBlocks;
}


// ************************************************************************* //
