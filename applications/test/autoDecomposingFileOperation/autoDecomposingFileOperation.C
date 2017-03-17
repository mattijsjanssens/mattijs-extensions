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

#include "autoDecomposingFileOperation.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvFieldDecomposer.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(autoDecomposingFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        autoDecomposingFileOperation,
        word
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::autoDecomposingFileOperation::
autoDecomposingFileOperation
(
    const bool verbose
)
:
    uncollatedFileOperation(false)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::autoDecomposingFileOperation::
~autoDecomposingFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::fileOperations::autoDecomposingFileOperation::baseMesh
(
    const Time& runTime
) const
{
    if (!baseMeshPtr_.valid())
    {
        InfoInFunction<< "Create time\n" << endl;
        baseRunTimePtr_.reset
        (
            new Time
            (
                runTime.controlDict(),
                runTime.rootPath(),
                runTime.globalCaseName(),
                runTime.system(),
                runTime.constant(),
                false                   // enableFunctionObjects
            )
        );

        InfoInFunction<< "Create mesh for time = "
            << runTime.timeName() << nl << endl;
        baseMeshPtr_.reset
        (
            new fvMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    baseRunTimePtr_().timeName(),
                    baseRunTimePtr_(),
                    IOobject::MUST_READ
                )
            )
        );
    }
    return baseMeshPtr_();
}


// const Foam::fvMesh&
// Foam::fileOperations::autoDecomposingFileOperation::procMesh
// (
//     const Time& runTime
// ) const
// {
//     if (!procMeshPtr_.valid())
//     {
//         InfoInFunction<< "Create mesh for time = "
//             << runTime.timeName() << nl << endl;
//         procMeshPtr_.reset
//         (
//             new fvMesh
//             (
//                 IOobject
//                 (
//                     fvMesh::defaultRegion,
//                     runTime.timeName(),
//                     runTime,
//                     IOobject::MUST_READ
//                 )
//             )
//         );
//     }
//     return procMeshPtr_();
// }


Foam::fileName Foam::fileOperations::autoDecomposingFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io
) const
{
    if (debug)
    {
        Pout<< "autoDecomposingFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }


    fileName objPath;
    if (io.name() == "p")
    {
        Pout<< "** searching global for p" << endl;

        fileName path = io.path();
        fileName objectPath = path/io.name();

        if (isFile(objectPath))
        {
            objPath = objectPath;
        }
        else if (io.time().processorCase())
        {
            fileName parentObjectPath =
                io.rootPath()/io.time().globalCaseName()
               /io.instance()/io.db().dbDir()/io.local()/io.name();

            if (isFile(parentObjectPath))
            {
                objPath = parentObjectPath;
            }
        }

        if (objPath.empty())
        {
            objPath = uncollatedFileOperation::filePath(true, io);
        }
    }
    else
    {
        objPath = uncollatedFileOperation::filePath(checkGlobal, io);
    }

    //if (debug)
    {
        Pout<< "autoDecomposingFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::autoPtr<Foam::ISstream>
Foam::fileOperations::autoDecomposingFileOperation::readStream
(
    regIOobject& io,
    const fileName& fName,
    const bool valid
) const
{
DebugVar(io.name());
DebugVar(fName);



    return uncollatedFileOperation::readStream(io, fName, valid);
}


bool Foam::fileOperations::autoDecomposingFileOperation::read
(
    regIOobject& io,
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& type
) const
{
    bool ok = true;

DebugVar(type);
DebugVar(io.name());


    if (type == volScalarField::typeName && Pstream::parRun())
    {
        Pout<< "** reading volScalarField " << io.name() << endl;

        // Set flag for e.g. codeStream
        const bool oldGlobal = io.globalObject();
        io.globalObject() = masterOnly;
        // If codeStream originates from dictionary which is
        // not IOdictionary we have a problem so use global
//        const bool oldFlag = regIOobject::masterOnlyReading;
//        regIOobject::masterOnlyReading = masterOnly;


        const fvMesh& procMesh = dynamic_cast<const fvMesh&>(io.db());
        // Force reading of the baseMesh
        const fvMesh& undecomposedMesh = baseMesh(io.time());

        // Force reading of the decomposition maps
        if (decomposerPtr_.valid())
        {
            InfoInFunction<< "Loading maps\n" << endl;
            cellAddressingPtr_.reset
            (
                new labelIOList
                (
                    IOobject
                    (
                        "cellProcAddressing",
                        procMesh.facesInstance(),
                        procMesh.meshSubDir,
                        procMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            );
            faceAddressingPtr_.reset
            (
                new labelIOList
                (
                    IOobject
                    (
                        "faceProcAddressing",
                        procMesh.facesInstance(),
                        procMesh.meshSubDir,
                        procMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            );
            boundaryAddressingPtr_.reset
            (
                new labelIOList
                (
                    IOobject
                    (
                        "boundaryProcAddressing",
                        procMesh.facesInstance(),
                        procMesh.meshSubDir,
                        procMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                )
            );
            InfoInFunction<< "Constructing decomposer\n" << endl;
            decomposerPtr_.reset
            (
                new fvFieldDecomposer
                (
                    undecomposedMesh,
                    procMesh,
                    faceAddressingPtr_(),
                    cellAddressingPtr_(),
                    boundaryAddressingPtr_()
                )
            );
        }

        // Read file
        ok = io.readData(io.readStream(typeName));
        io.close();

        // Restore flags
        io.globalObject() = oldGlobal;
//        regIOobject::masterOnlyReading = oldFlag;
    }
    else
    {
        ok = uncollatedFileOperation::read
        (
            io,
            masterOnly,
            format,
            type
        );
    }

    return ok;
}


// ************************************************************************* //
