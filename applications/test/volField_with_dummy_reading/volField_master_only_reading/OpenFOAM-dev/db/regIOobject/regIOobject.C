/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "regIOobject.H"
#include "Time.H"
#include "polyMesh.H"
#include "registerSwitch.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regIOobject, 0);

    template<>
    const char* NamedEnum
    <
        regIOobject::fileCheckTypes,
        4
    >::names[] =
    {
        "timeStamp",
        "timeStampMaster",
        "inotify",
        "inotifyMaster"
    };
}

float Foam::regIOobject::fileModificationSkew
(
    Foam::debug::floatOptimisationSwitch("fileModificationSkew", 30)
);
registerOptSwitch
(
    "fileModificationSkew",
    float,
    Foam::regIOobject::fileModificationSkew
);


const Foam::NamedEnum<Foam::regIOobject::fileCheckTypes, 4>
    Foam::regIOobject::fileCheckTypesNames;

Foam::regIOobject::fileCheckTypes Foam::regIOobject::fileModificationChecking
(
    fileCheckTypesNames.read
    (
        debug::optimisationSwitches().lookup
        (
            "fileModificationChecking"
        )
    )
);

namespace Foam
{
    // Register re-reader
    class addfileModificationCheckingToOpt
    :
        public ::Foam::simpleRegIOobject
    {
    public:

        addfileModificationCheckingToOpt(const char* name)
        :
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
        {}

        virtual ~addfileModificationCheckingToOpt()
        {}

        virtual void readData(Foam::Istream& is)
        {
            regIOobject::fileModificationChecking =
                regIOobject::fileCheckTypesNames.read(is);
        }

        virtual void writeData(Foam::Ostream& os) const
        {
            os <<  regIOobject::fileCheckTypesNames
                [regIOobject::fileModificationChecking];
        }
    };

    addfileModificationCheckingToOpt addfileModificationCheckingToOpt_
    (
        "fileModificationChecking"
    );
}


bool Foam::regIOobject::masterOnlyReading = false;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regIOobject::regIOobject(const IOobject& io, const bool isTime)
:
    IOobject(io),
    registered_(false),
    ownedByRegistry_(false),
    watchIndex_(-1),
    eventNo_                // Do not get event for top level Time database
    (
        isTime
      ? 0
      : db().getEvent()
    )
{
    // Register with objectRegistry if requested
    if (registerObject())
    {
        checkIn();
    }
}


Foam::regIOobject::regIOobject(const regIOobject& rio)
:
    IOobject(rio),
    registered_(false),
    ownedByRegistry_(false),
    watchIndex_(rio.watchIndex_),
    eventNo_(db().getEvent())
{
    // Do not register copy with objectRegistry
}


Foam::regIOobject::regIOobject(const regIOobject& rio, bool registerCopy)
:
    IOobject(rio),
    registered_(false),
    ownedByRegistry_(false),
    watchIndex_(-1),
    eventNo_(db().getEvent())
{
    if (registerCopy && rio.registered_)
    {
        const_cast<regIOobject&>(rio).checkOut();
        checkIn();
    }
}


Foam::regIOobject::regIOobject
(
    const word& newName,
    const regIOobject& rio,
    bool registerCopy
)
:
    IOobject(newName, rio.instance(), rio.local(), rio.db()),
    registered_(false),
    ownedByRegistry_(false),
    watchIndex_(-1),
    eventNo_(db().getEvent())
{
    if (registerCopy)
    {
        checkIn();
    }
}


Foam::regIOobject::regIOobject
(
    const IOobject& io,
    const regIOobject& rio
)
:
    IOobject(io),
    registered_(false),
    ownedByRegistry_(false),
    watchIndex_(-1),
    eventNo_(db().getEvent())
{
    if (registerObject())
    {
        checkIn();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regIOobject::~regIOobject()
{
    if (objectRegistry::debug)
    {
        Info<< "Destroying regIOobject called " << name()
            << " of type " << type()
            << " in directory " << path()
            << endl;
    }

    // Check out of objectRegistry if not owned by the registry

    if (!ownedByRegistry_)
    {
        checkOut();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regIOobject::checkIn()
{
    if (!registered_)
    {
        // multiple checkin of same object is disallowed - this would mess up
        // any mapping
        registered_ = db().checkIn(*this);

        if
        (
            registered_
         && readOpt() == MUST_READ_IF_MODIFIED
         && time().runTimeModifiable()
        )
        {
            if (watchIndex_ != -1)
            {
                FatalErrorInFunction
                    << "Object " << objectPath()
                    << " already watched with index " << watchIndex_
                    << abort(FatalError);
            }

            fileName f = fileHandler().filePath(*this);
            if (!f.size())
            {
                // We don't have this file but would like to re-read it.
                // Possibly if master-only reading mode.
                f = objectPath();
            }
            watchIndex_ = time().addWatch(f);
        }

        // check-in on defaultRegion is allowed to fail, since subsetted meshes
        // are created with the same name as their originating mesh
        if (!registered_ && debug && name() != polyMesh::defaultRegion)
        {
            if (debug == 2)
            {
                // for ease of finding where attempted duplicate check-in
                // originated
                FatalErrorInFunction
                    << "failed to register object " << objectPath()
                    << " the name already exists in the objectRegistry" << endl
                    << "Contents:" << db().sortedToc()
                    << abort(FatalError);
            }
            else
            {
                WarningInFunction
                    << "failed to register object " << objectPath()
                    << " the name already exists in the objectRegistry"
                    << endl;
            }
        }
    }

    return registered_;
}


bool Foam::regIOobject::checkOut()
{
    if (registered_)
    {
        registered_ = false;

        if (watchIndex_ != -1)
        {
            time().removeWatch(watchIndex_);
            watchIndex_ = -1;
        }
        return db().checkOut(*this);
    }

    return false;
}


bool Foam::regIOobject::upToDate(const regIOobject& a) const
{
    if (a.eventNo() >= eventNo_)
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::regIOobject::upToDate
(
    const regIOobject& a,
    const regIOobject& b
) const
{
    if
    (
        a.eventNo() >= eventNo_
     || b.eventNo() >= eventNo_
    )
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::regIOobject::upToDate
(
    const regIOobject& a,
    const regIOobject& b,
    const regIOobject& c
) const
{
    if
    (
        a.eventNo() >= eventNo_
     || b.eventNo() >= eventNo_
     || c.eventNo() >= eventNo_
    )
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::regIOobject::upToDate
(
    const regIOobject& a,
    const regIOobject& b,
    const regIOobject& c,
    const regIOobject& d
) const
{
    if
    (
        a.eventNo() >= eventNo_
     || b.eventNo() >= eventNo_
     || c.eventNo() >= eventNo_
     || d.eventNo() >= eventNo_
    )
    {
        return false;
    }
    else
    {
        return true;
    }
}


void Foam::regIOobject::setUpToDate()
{
    eventNo_ = db().getEvent();
}


void Foam::regIOobject::rename(const word& newName)
{
    // Check out of objectRegistry
    checkOut();

    IOobject::rename(newName);

    if (registerObject())
    {
        // Re-register object with objectRegistry
        checkIn();
    }
}


void Foam::regIOobject::operator=(const IOobject& io)
{
    isPtr_.clear();

    // Check out of objectRegistry
    checkOut();

    IOobject::operator=(io);

    if (registerObject())
    {
        // Re-register object with objectRegistry
        checkIn();
    }
}


// ************************************************************************* //
