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

Class
    Foam::fileOperations::autoParallelFileOperation

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_autoParallelFileOperation_H
#define fileOperations_autoParallelFileOperation_H

#include "uncollatedFileOperation.H"
#include "IOstreams.H"
#include "collatedFileOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class unallocatedFvMesh;
class mapDistributePolyMesh;
class parUnallocatedFvFieldReconstructor;

namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                 Class autoParallelFileOperation Declaration
\*---------------------------------------------------------------------------*/

class autoParallelFileOperation
:
    public uncollatedFileOperation
{
public:

    class storeFileHandler
    {
        autoPtr<fileOperation>& basicHandler_;

        public:

            storeFileHandler(autoPtr<fileOperation>& basicHandler)
            :
                basicHandler_(basicHandler)
            {
                if
                (
                    basicHandler().type()
                 != fileOperation::fileHandlerPtr_().type()
                )
                {
                    fileOperation* p = fileOperation::fileHandlerPtr_.ptr();

                    if (autoParallelFileOperation::debug)
                    {
                        Pout<< "Saving fileHandler:" << p->type()
                            << " installing fileHandler:"
                            << basicHandler().type() << endl;
                    }

                    (void)fileHandler(basicHandler_);

                    basicHandler_.set(p);
                }
            }

            ~storeFileHandler()
            {
                if
                (
                    basicHandler_().type()
                 != fileOperation::fileHandlerPtr_().type()
                )
                {
                    if (autoParallelFileOperation::debug)
                    {
                        Pout<< "Restoring saved fileHandler:"
                            << basicHandler_().type() << endl;
                    }

                    fileOperation* p = basicHandler_.ptr();
                    basicHandler_.set(fileOperation::fileHandlerPtr_.ptr());
                    autoPtr<fileOperation> restore(p);
                    (void)fileHandler(restore);
                }
            }
    };



    // Parallel on-the-fly reconstruction

        //- Current file handler
        mutable autoPtr<fileOperation> basicFileHandler_;

        //- Proc faces instance (where procAddressing gets read from)
        mutable word procFacesInstance_;

        //- Undecomposed database
        mutable autoPtr<Time> baseRunTimePtr_;

        //- Decomposition
        mutable autoPtr<mapDistributePolyMesh> distMapPtr_;

        //- Undecomposed case (basic info only; no addressing)
        mutable autoPtr<unallocatedFvMesh> baseMeshPtr_;

        //- Processor mesh (basic info only; no addressing)
        mutable autoPtr<unallocatedFvMesh> meshPtr_;

        //- Reconstruction engine
        mutable autoPtr<parUnallocatedFvFieldReconstructor> reconstructorPtr_;


    // Private Member Functions

        //- Return the base Time (demand loaded)
        const Time& baseRunTime(const Time&) const;

        //- Return the base mesh (demand loaded)
        const unallocatedFvMesh& baseMesh(const Time&) const;

        //- Return the current (decomposed) mesh (demand loaded)
        const unallocatedFvMesh& mesh(const Time&) const;


public:

        //- Runtime type information
        TypeName("autoParallel");


    // Constructors

        //- Construct null
        autoParallelFileOperation(const bool verbose);


    //- Destructor
    virtual ~autoParallelFileOperation();


    // Member Functions

        // fileOperation functinality

            //- Search for a file or directory. Use IOobject version in
            //  preference
            //using fileOperation::filePath;
            virtual fileName filePath(const fileName&) const;

            //- Search for an object. checkGlobal : also check undecomposed case
            //  This differs from uncollated in searching for parent
            //  also for files in time directories
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&,
                const word& typeName
            ) const;

            //- Get sorted list of times
            virtual instantList findTimes(const fileName&, const word&) const;

            //- Callback for time change
            virtual void setTime(const Time&) const;

            //- Search directory for objects. Used in IOobjectList.
            virtual fileNameList readObjects
            (
                const objectRegistry& db,
                const fileName& instance,
                const fileName& local,
                word& newInstance
            ) const;

            //- Top-level read
            virtual bool read
            (
                regIOobject&,
                const bool masterOnly,
                const IOstream::streamFormat format,
                const word& typeName
            ) const;

            //- Writes a regIOobject
            virtual bool writeObject
            (
                const regIOobject&,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED,
                const bool valid = true
            ) const;
};


/*---------------------------------------------------------------------------*\
             Class autoParallelFileOperationInitialise Declaration
\*---------------------------------------------------------------------------*/

// class autoParallelFileOperationInitialise
// :
//     public collatedFileOperationInitialise
// {
// public:
//
//     // Constructors
//
//         //- Construct from components
//         autoParallelFileOperationInitialise(int& argc, char**& argv);
//
//
//     //- Destructor
//     virtual ~autoParallelFileOperationInitialise()
//     {}
// };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
