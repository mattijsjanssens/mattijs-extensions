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

#include "fileOperation.H"
#include "localFileOperation.H"
#include "regIOobject.H"
#include "argList.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    autoPtr<fileOperation> fileOperation::fileHandlerPtr_;

    defineTypeNameAndDebug(fileOperation, 0);
    defineRunTimeSelectionTable(fileOperation, word);

    class addArgsOptions
    {
        public:
        addArgsOptions()
        {
            argList::addOption
            (
                "fileHandler",
                "handler",
                "override the fileHandler"
            );
        }
    };

    addArgsOptions intObj;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperation::fileOperation()
{}


Foam::autoPtr<Foam::fileOperation> Foam::fileOperation::New(const word& type)
{
    if (debug)
    {
        InfoInFunction << "Constructing fileOperation" << endl;
    }

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(type);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fileOperation type "
            << type << nl << nl
            << "Valid fileOperation types are" << endl
            << wordConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<fileOperation>(cstrIter()());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperation::~fileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperation::writeObject
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


const Foam::fileOperation& Foam::fileHandler()
{
    if (!fileOperation::fileHandlerPtr_.valid())
    {
        word handler(getEnv("FOAM_FILEHANDLER"));
        if (!handler.size())
        {
            handler = fileOperations::localFileOperation::typeName;
        }

        cout<< "fileHandler() : Inserting fileOperation of type "
            << handler << std::endl;
        fileOperation::fileHandlerPtr_ = fileOperation::New(handler);
    }
    return fileOperation::fileHandlerPtr_();
}


void Foam::fileHandler(autoPtr<fileOperation>& newHandlerPtr)
{
    if (fileOperation::fileHandlerPtr_.valid())
    {
        cout<< "fileHandler() : Deleting fileOperation of type "
            << fileOperation::fileHandlerPtr_().type() << std::endl;
    }
    fileOperation::fileHandlerPtr_.clear();

    if (newHandlerPtr.valid())
    {
        cout<< "fileHandler() : Inserting fileOperation of type "
            << newHandlerPtr().type() << std::endl;
        fileOperation::fileHandlerPtr_ = newHandlerPtr;
    }
}


// ************************************************************************* //
