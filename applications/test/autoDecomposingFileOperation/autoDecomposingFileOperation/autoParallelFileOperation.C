/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "autoParallelFileOperation.H"
#include "Time.H"
#include "fvMesh.H"
//#include "fvFieldDecomposer.H"
#include "addToRunTimeSelectionTable.H"
//#include "processorMeshes.H"
//#include "fvFieldReconstructor.H"

#include "mapDistributePolyMesh.H"
#include "unallocatedFvMesh.H"
#include "unallocatedFvMeshTools.H"
#include "parUnallocatedFvFieldReconstructor.H"
#include "uVolFields.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(autoParallelFileOperation, 0);
    addRemovableToRunTimeSelectionTable
    (
        fileOperation,
        autoParallelFileOperation,
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
            //    new autoParallelFileOperation(true)
            //);
            //Foam::fileHandler(handler);
        }

        ~installFileOp()
        {
            if (fileHandler().type() == autoParallelFileOperation::typeName)
            {
                autoPtr<fileOperation> handler(nullptr);
                Foam::fileHandler(handler);
            }
        }
    };
    installFileOp installAutoParallel_;
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Time& Foam::fileOperations::autoParallelFileOperation::baseRunTime
(
    const Time& runTime
) const
{
    if (baseRunTimePtr_.valid())
    {
        Info<< "Creating base time\n" << endl;
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
    }
    baseRunTimePtr_().setTime(runTime);
    return baseRunTimePtr_();
}


const Foam::unallocatedFvMesh&
Foam::fileOperations::autoParallelFileOperation::baseMesh
(
    const Time& runTime
) const
{
    if (!baseMeshPtr_.valid())
    {
        Info<< "Reading decomposition addressing" << nl << endl;
        distMapPtr_ = unallocatedFvMeshTools::readReconstructMap
        (
            IOobject
            (
                "dummy",
                runTime.findInstance(fvMesh::meshSubDir, "faces"),
                fvMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        Info<< "Creating base mesh for time = "
            << runTime.timeName() << nl << endl;
        baseMeshPtr_ = unallocatedFvMeshTools::newMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,      // name of mesh
                baseRunTime(runTime).timeName(),
                baseRunTime(runTime),
                IOobject::MUST_READ
            ),
            distMapPtr_().cellMap().constructSize()
        );
    }
    return baseMeshPtr_();
}


// const Foam::parUnallocatedFvFieldReconstructor&
// Foam::fileOperations::autoParallelFileOperation::reconstructor
// (
//     const unallocatedFvMesh& mesh
// ) const
// {
//     if (!reconstructorPtr_.valid())
//     {
//         Info<< "Creating reconstructor" << nl << endl;
//         reconstructorPtr_ = new parUnallocatedFvFieldReconstructor
//         (
//             baseMesh(),
//             mesh,
//             distMapPtr_()
//         );
//     }
//     return reconstructorPtr_();
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::autoParallelFileOperation::
autoParallelFileOperation
(
    const bool verbose
)
:
    uncollatedFileOperation(false)
{
    if (verbose)
    {
        Info<< "I/O    : " << typeName << nl
            << "         (reconstructs on-the-fly in parallel operation)"
            << endl;
    }

    if (regIOobject::fileModificationChecking == regIOobject::timeStampMaster)
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
        regIOobject::fileModificationChecking = regIOobject::timeStamp;
    }
    else if
    (
        regIOobject::fileModificationChecking
     == regIOobject::inotifyMaster
    )
    {
        if (verbose)
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify"
                << endl;
        }
        regIOobject::fileModificationChecking = regIOobject::inotify;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::autoParallelFileOperation::
~autoParallelFileOperation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Foam::fileName Foam::fileOperations::autoParallelFileOperation::filePath
// (
//     const bool checkGlobal,
//     const IOobject& io,
//     const word& typeName
// ) const
// {
//     if (debug)
//     {
//         Pout<< "autoParallelFileOperation::filePath :"
//             << " objectPath:" << io.objectPath()
//             << " checkGlobal:" << checkGlobal << endl;
//     }
// 
//     // Try uncollated searching
//     fileName objPath = uncollatedFileOperation::filePath
//     (
//         checkGlobal,
//         io,
//         typeName
//     );
// 
//     // If not found and parallel check parent
//     if (objPath.empty() && checkGlobal && io.time().processorCase())
//     {
//         fileName parentObjectPath =
//             io.rootPath()/io.time().globalCaseName()
//            /io.instance()/io.db().dbDir()/io.local()/io.name();
// 
//         if (isFile(parentObjectPath))
//         {
//             objPath = parentObjectPath;
//         }
//     }
// 
//     if (debug)
//     {
//         Pout<< "autoParallelFileOperation::filePath :"
//             << " Returning from file searching:" << endl
//             << "    objectPath:" << io.objectPath() << endl
//             << "    filePath  :" << objPath << endl << endl;
//     }
//     return objPath;
// }
// 
// 
// bool Foam::fileOperations::autoParallelFileOperation::read
// (
//     regIOobject& io,
//     const bool masterOnly,
//     const IOstream::streamFormat format,
//     const word& type
// ) const
// {
//     bool ok = true;
// 
//     if
//     (
//         Pstream::parRun()
//      && (
//             type == volScalarField::typeName
//          || type == volVectorField::typeName
//          || type == volSphericalTensorField::typeName
//          || type == volSymmTensorField::typeName
//          || type == volTensorField::typeName
//         )
//     )
//     {
//         // Set flag for e.g. codeStream
//         const bool oldGlobal = io.globalObject();
//         io.globalObject() = masterOnly;
//         // If codeStream originates from dictionary which is
//         // not IOdictionary we have a problem so use global
//         //const bool oldFlag = regIOobject::masterOnlyReading;
//         //regIOobject::masterOnlyReading = masterOnly;
// 
// 
//         // Find file, check in parent directory
//         fileName objPath = filePath(true, io, type);
// 
//         // Check if the file comes from the parent path
//         fileName parentObjectPath =
//             io.rootPath()/io.time().globalCaseName()
//            /io.instance()/io.db().dbDir()/io.local()/io.name();
// 
//         if (debug)
//         {
//             Pout<< "io.objectPath   :" << io.objectPath() << endl;
//             Pout<< "filePath        :" << objPath << endl;
//             Pout<< "parentObjectPath:" << parentObjectPath << endl;
//         }
// 
//         if (io.objectPath() != objPath && objPath == parentObjectPath)
//         {
//             // Force reading of the baseMesh
//             const fvMesh& undecomposedMesh = baseMesh(io.time());
// 
//             IOobject parentIO
//             (
//                 io.name(),
//                 io.instance(),
//                 io.local(),
//                 undecomposedMesh,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE,
//                 false
//             );
// 
//             OStringStream os;
// 
//             decomposeAndWrite<volScalarField>(io, parentIO, type, os);
//             decomposeAndWrite<volVectorField>(io, parentIO, type, os);
//             decomposeAndWrite<volSphericalTensorField>(io, parentIO, type, os);
//             decomposeAndWrite<volSymmTensorField>(io, parentIO, type, os);
//             decomposeAndWrite<volTensorField>(io, parentIO, type, os);
// 
//             IStringStream is(os.str());
// 
//             // Read field from stream
//             ok = io.readData(is);
//             io.close();
//         }
//         else
//         {
//             ok = io.readData(io.readStream(type));
//             io.close();
//         }
// 
//         // Restore flags
//         io.globalObject() = oldGlobal;
//         //regIOobject::masterOnlyReading = oldFlag;
//     }
//     else
//     {
//         ok = uncollatedFileOperation::read
//         (
//             io,
//             masterOnly,
//             format,
//             type
//         );
//     }
// 
//     return ok;
// }


bool Foam::fileOperations::autoParallelFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    bool ok = true;

    if
    (
        Pstream::parRun()
     && (
            io.type() == volScalarField::typeName
         //|| io.type() == volVectorField::typeName
         //|| io.type() == volSphericalTensorField::typeName
         //|| io.type() == volSymmTensorField::typeName
         //|| io.type() == volTensorField::typeName
        )
    )
    {
        if (!valid)
        {
            return ok;
        }

Pout<< "** reconstructing:" << io.objectPath() << endl;
        const volScalarField& fld = dynamic_cast<const uVolScalarField&>(io);
Pout<< "** reconstructing:" << fld.name() << endl;
        return reconstructAndWrite(fld, fmt, ver, cmp);
    }
    else
    {
Pout<< "** non-par:" << io.objectPath() << endl;

        ok = uncollatedFileOperation::writeObject
        (
            io,
            fmt,
            ver,
            cmp,
            valid
        );
    }
    return ok;
}


// ************************************************************************* //
