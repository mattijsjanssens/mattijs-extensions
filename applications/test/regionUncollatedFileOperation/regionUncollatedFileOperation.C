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

#include "regionUncollatedFileOperation.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"
#include "IOdictionary.H"
#include "volFields.H"
//#include "fieldDictionary.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(regionUncollatedFileOperation, 0);
    addRemovableToRunTimeSelectionTable
    (
        fileOperation,
        regionUncollatedFileOperation,
        word
    );

    class installFileOp
    {
    public:

        installFileOp()
        {
            // Install autoDecomposing as fileHandler
            //autoPtr<fileOperation> handler
            //(
            //    new regionUncollatedFileOperation(true)
            //);
            //Foam::fileHandler(handler);
        }

        ~installFileOp()
        {
            if (fileHandler().type() == regionUncollatedFileOperation::typeName)
            {
                autoPtr<fileOperation> handler(nullptr);
                Foam::fileHandler(handler);
            }
        }
    };
    installFileOp installFileOp_;
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::regionUncollatedFileOperation::
regionUncollatedFileOperation(const bool verbose)
:
    uncollatedFileOperation(false)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::regionUncollatedFileOperation::
~regionUncollatedFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileNameList
Foam::fileOperations::regionUncollatedFileOperation::readObjects
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    word& newInstance
) const
{
    if (debug)
    {
        Pout<< "regionUncollatedFileOperation::readObjects :"
            << " db:" << db.objectPath()
            << " instance:" << instance << endl;
    }

    //- Use non-time searching version
    fileNameList objectNames
    (
        uncollatedFileOperation::readObjects
        (
            db,
            instance,
            local,
            newInstance
        )
    );

    if (objectNames.empty())
    {
        // Did not find files.  Try directory up. Use either parent of local
        // or parent of the database.
        if (!local.empty())
        {
            return uncollatedFileOperation::readObjects
            (
                db,
                instance,
                local.path(),
                newInstance
            );
        }
        else if (&db != dynamic_cast<const objectRegistry*>(&db.time()))
        {
            return uncollatedFileOperation::readObjects
            (
                db.parent(),
                instance,
                local.path(),
                newInstance
            );
        }
    }

    if (debug)
    {
        Pout<< "regionUncollatedFileOperation::readObjects :"
            << " newInstance:" << newInstance
            << " objectNames:" << objectNames << endl;
    }

    return objectNames;
}


Foam::fileName Foam::fileOperations::regionUncollatedFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& type
) const
{
    if (debug)
    {
        Pout<< "regionUncollatedFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal
            << " type:" << type << endl;
    }

    // Try uncollated searching
    fileName objPath = uncollatedFileOperation::filePath(checkGlobal, io, type);

    if (objPath.empty())
    {
        // Try directory up. Use either parent of local or parent of the
        // database.
        if (!io.local().empty())
        {
            IOobject parentIO(io);
            const_cast<fileName&>(parentIO.local()) = parentIO.local().path();
            if (debug)
            {
                Pout<< "Trying parent-of-local:" << io.local()
                    << " objPath:" << parentIO.objectPath() << endl;
            }
            objPath = uncollatedFileOperation::filePath
            (
                checkGlobal,
                parentIO,
                type
            );
        }
        else if
        (
           &io.db()
         != dynamic_cast<const objectRegistry*>(&io.time())
        )
        {
            // Use the parent of the database
            IOobject parentIO(io, io.db().parent());
            if (debug)
            {
                Pout<< "Trying parent-of-database objPath:"
                    << parentIO.objectPath() << endl;
            }
            objPath = uncollatedFileOperation::filePath
            (
                checkGlobal,
                parentIO,
                type
            );
        }
    }

    if (debug)
    {
        Pout<< "regionUncollatedFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


bool Foam::fileOperations::regionUncollatedFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& type
) const
{
    if (debug)
    {
        Pout<< "regionUncollatedFileOperation::read :"
            << " objectPath:" << io.objectPath()
            << " type:" << type << endl;
    }

    static const regExp volFields("vol.*Field");
    static const regExp uniformDimensionedFields("uniformDimensioned.*Field");

    if
    (
        type != IOdictionary::typeName
     && !volFields.match(type)
     && !uniformDimensionedFields.match(type)
    )
    {
        if (debug)
        {
            Pout<< "Normal uncollatedFileOperation::read" << endl;
        }
        return uncollatedFileOperation::read
        (
            io,
            masterOnly,
            format,
            type
        );
    }

    bool checkGlobal = (type == IOdictionary::typeName);


    // Try to locate the dictionary
    fileName objPath = uncollatedFileOperation::filePath
    (
        checkGlobal,
        io,
        type
    );

    bool ok = false;
    if (!objPath.empty())
    {
        if (debug)
        {
            Pout<< "Reverting to uncollatedFileOperation::read" << endl;
        }
        ok = uncollatedFileOperation::read
        (
            io,
            masterOnly,
            format,
            type
        );
    }
    else
    {
        // Did not find file.  Try directory up. Use either parent of local
        // or parent of the database.

        word regionName;
        autoPtr<IOobject> parentIOPtr;

        if (!io.local().empty())
        {
            // Modifiy the local field
            parentIOPtr.reset(new IOobject(io));
            IOobject& parentIO = parentIOPtr();
            regionName = parentIO.local().name();
            const_cast<fileName&>(parentIO.local()) = parentIO.local().path();

            if (debug)
            {
                Pout<< "Trying parent-of-local:" << io.local()
                    << " region:" << regionName
                    << " objPath:" << parentIO.objectPath() << endl;
            }
        }
        else if
        (
           &io.db()
         != dynamic_cast<const objectRegistry*>(&io.time())
        )
        {
            // Use the parent of the database
            parentIOPtr.reset(new IOobject(io, io.db().parent()));
            IOobject& parentIO = parentIOPtr();
            regionName = io.db().name();

            if (debug)
            {
                Pout<< "Trying parent-of-database objPath:"
                    << " region:" << regionName
                    << " objPath:" << parentIO.objectPath() << endl;
            }
        }


        //if (data.found(regionName))
        if (!regionName.empty())
        {
            // Load parent dictionary
            dictionary data;
            if (checkGlobal)
            {
                IOdictionary parentDict(parentIOPtr(), type);
                data.transfer(parentDict);
            }
            else
            {
                localIOdictionary parentDict(parentIOPtr(), type);
                data.transfer(parentDict);
            }

            OStringStream os;
            os << data.subDict(regionName);
            IStringStream is(os.str());
            ok = io.readData(is);
        }
    }

    return ok;
}


// ************************************************************************* //
