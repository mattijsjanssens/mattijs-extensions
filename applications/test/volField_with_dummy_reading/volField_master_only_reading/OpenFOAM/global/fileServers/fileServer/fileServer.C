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

#include "fileServer.H"
//#include "masterFileServer.H"
#include "localFileServer.H"
#include "regIOobject.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    autoPtr<fileServer> fileServer::serverPtr_;

    defineTypeNameAndDebug(fileServer, 0);
    defineRunTimeSelectionTable(fileServer, word);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileServer::fileServer()
{}


Foam::autoPtr<Foam::fileServer> Foam::fileServer::New(const word& type)
{
    if (debug)
    {
        InfoInFunction << "Constructing fileServer" << endl;
    }

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(type);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fileServer type "
            << type << nl << nl
            << "Valid fileServer types are" << endl
            << wordConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<fileServer>(cstrIter()());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileServer::~fileServer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileServer::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    mkDir(io.path());

    fileName pathName(io.objectPath());

    autoPtr<Ostream> osPtr
    (
        NewOFstream
        (
            pathName,
            fmt,
            ver,
            cmp
        )
    );

    if (!osPtr.valid())
    {
        return false;
    }

    Ostream& os = osPtr();

    // If any of these fail, return (leave error handling to Ostream class)
    if (!os.good())
    {
        return false;
    }

    if (!io.writeHeader(os))
    {
        return false;
    }

    // Write the data to the Ostream
    if (!io.writeData(os))
    {
        return false;
    }

    IOobject::writeEndDivider(os);

    return true;
}


const Foam::fileServer& Foam::server()
{
    if (!Foam::fileServer::serverPtr_.valid())
    {
        cout<< "Foam::server() : Allocating fileServer" << std::endl;

        Foam::fileServer::serverPtr_.reset(new fileServers::localFileServer());
    }
    return Foam::fileServer::serverPtr_();
}


const Foam::fileServer& Foam::server(autoPtr<fileServer>& newServerPtr)
{
    if (Foam::fileServer::serverPtr_.valid())
    {
        cout<< "Foam::server() : Deleting fileServer of type "
            << Foam::fileServer::serverPtr_().type() << std::endl;
    }
    Foam::fileServer::serverPtr_.clear();

    if (newServerPtr.valid())
    {
        cout<< "Foam::server() : Inserting fileServer of type "
            << newServerPtr().type() << std::endl;
        Foam::fileServer::serverPtr_ = newServerPtr;
    }
    return Foam::fileServer::serverPtr_();
}


// ************************************************************************* //
