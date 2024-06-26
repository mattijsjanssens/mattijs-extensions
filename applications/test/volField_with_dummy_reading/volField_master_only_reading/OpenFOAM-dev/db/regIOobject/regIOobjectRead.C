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
#include "IFstream.H"
#include "Time.H"
#include "Pstream.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::regIOobject::readStream()
{
    if (IFstream::debug)
    {
        //InfoInFunction
        Pout<< "regIOobject::readStream() : "
            << "Reading object " << name()
            << " from file " << objectPath()
            << endl;
    }

    if (readOpt() == NO_READ)
    {
        FatalErrorInFunction
            << "NO_READ specified for read-constructor of object " << name()
            << " of class " << headerClassName()
            << abort(FatalError);
    }

    // Construct object stream and read header if not already constructed
    if (!isPtr_.valid())
    {
        fileName objPath;
        if (watchIndex_ != -1)
        {
            // File is being watched. Read exact file that is being watched.
            objPath = time().getFile(watchIndex_);
        }
        else
        {
            // Search intelligently for file
            //objPath = filePath();
            objPath = fileHandler().filePath(*this);

            if (!objPath.size())
            {
                FatalIOError
                (
                    "regIOobject::readStream()",
                    __FILE__,
                    __LINE__,
                    objectPath(),
                    0
                )   << "cannot find file"
                    << exit(FatalIOError);
            }
        }

        isPtr_ = fileHandler().readStream(*this, objPath);
    }

    // Mark as uptodate if read successfully
    if (watchIndex_ != -1)
    {
        time().setUnmodified(watchIndex_);
    }

    return isPtr_();
}


Foam::Istream& Foam::regIOobject::readStream(const word& expectName)
{
    if (IFstream::debug)
    {
        InfoInFunction
            << "Reading object " << name()
            << " from file " << objectPath()
            << endl;
    }

    // Construct IFstream if not already constructed
    if (!isPtr_.valid())
    {
        readStream();

        // Check the className of the regIOobject
        // dictionary is an allowable name in case the actual class
        // instantiated is a dictionary
        if
        (
            expectName.size()
         && headerClassName() != expectName
         && headerClassName() != "dictionary"
        )
        {
            FatalIOErrorInFunction(isPtr_())
                << "unexpected class name " << headerClassName()
                << " expected " << expectName << endl
                << "    while reading object " << name()
                << exit(FatalIOError);
        }
    }

    return isPtr_();
}


void Foam::regIOobject::close()
{
    if (IFstream::debug)
    {
        InfoInFunction
            << "Finished reading " << filePath()
            << endl;
    }

    isPtr_.clear();
}


bool Foam::regIOobject::readData(Istream&)
{
    return false;
}


bool Foam::regIOobject::read()
{
    // Note: cannot do anything in readStream itself since this is used by
    // e.g. GeometricField.

    bool masterOnly =
        regIOobject::fileModificationChecking == timeStampMaster
     || regIOobject::fileModificationChecking == inotifyMaster;

    bool ok = true;
    if (Pstream::master() || !masterOnly)
    {
        if (IFstream::debug)
        {
            Pout<< "regIOobject::read() : "
                << "reading object " << name()
                << " from file " << endl;
        }

        // Set flag for e.g. codeStream
        bool oldFlag = regIOobject::masterOnlyReading;
        regIOobject::masterOnlyReading = masterOnly;

        // Read file
        ok = readData(readStream(type()));
        close();

        regIOobject::masterOnlyReading = oldFlag;
    }

    if (masterOnly && Pstream::parRun())
    {
        // Scatter master data using communication scheme

        const List<Pstream::commsStruct>& comms =
        (
            (Pstream::nProcs() < Pstream::nProcsSimpleSum)
          ? Pstream::linearCommunication()
          : Pstream::treeCommunication()
        );

        // Master reads headerclassname from file. Make sure this gets
        // transfered as well as contents.
        Pstream::scatter
        (
            comms,
            const_cast<word&>(headerClassName()),
            Pstream::msgType(),
            Pstream::worldComm
        );
        Pstream::scatter(comms, note(), Pstream::msgType(), Pstream::worldComm);


        // Get my communication order
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            if (IFstream::debug)
            {
                Pout<< "regIOobject::read() : "
                    << "reading object " << name()
                    << " from processor " << myComm.above()
                    << endl;
            }

            // Note: use ASCII for now - binary IO of dictionaries is
            // not currently supported
            IPstream fromAbove
            (
                Pstream::scheduled,
                myComm.above(),
                0,
                Pstream::msgType(),
                Pstream::worldComm,
                IOstream::ASCII
            );
            ok = readData(fromAbove);
        }

        // Send to my downstairs neighbours. Note reverse order not
        // nessecary here - just for consistency reasons.
        forAllReverse(myComm.below(), belowI)
        {
            OPstream toBelow
            (
                Pstream::scheduled,
                myComm.below()[belowI],
                0,
                Pstream::msgType(),
                Pstream::worldComm,
                IOstream::ASCII
            );
            writeData(toBelow);
        }
    }
    return ok;
}


bool Foam::regIOobject::modified() const
{
    if (watchIndex_ != -1)
    {
        return time().getState(watchIndex_) != fileMonitor::UNMODIFIED;
    }
    else
    {
        return false;
    }
}


bool Foam::regIOobject::readIfModified()
{
    if (watchIndex_ != -1)
    {
        if (modified())
        {
            const fileName& fName = time().getFile(watchIndex_);
            InfoInFunction
                << nl
                << "    Re-reading object " << name()
                << " from file " << fName << endl;
            return read();
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
