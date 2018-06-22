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
    // Parallel on-the-fly reconstruction

        //- Undecomposed database
        mutable autoPtr<Time> baseRunTimePtr_;

        //- Decomposition
        mutable autoPtr<mapDistributePolyMesh> distMapPtr_;

        //- Undecomposed case (basic info only; no addressing)
        mutable autoPtr<unallocatedFvMesh> baseMeshPtr_;

        //- Reconstruction engine
        mutable autoPtr<parUnallocatedFvFieldReconstructor> reconstructorPtr_;


    // Private Member Functions

        //- Return the base Time (demand loaded)
        const Time& baseRunTime(const Time&) const;

        //- Return the base mesh (demand loaded)
        const unallocatedFvMesh& baseMesh(const Time&) const;

        //- Reconstruct and write
        template<class GeoField>
        bool reconstructAndWrite
        (
            const unallocatedFvMesh& uMesh,
            const GeoField& procFld,
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;

        //- Reconstruct and write
        template<class GeoField>
        bool reconstructAndWrite
        (
            const regIOobject& io,
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;

        //- Reconstruct and write a surface field
        template<class GeoField>
        bool reconstructAndWrite2
        (
            const regIOobject& io,
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;


public:

        //- Runtime type information
        TypeName("autoParallel");


    // Constructors

        //- Construct null
        autoParallelFileOperation(const bool verbose);


    //- Destructor
    virtual ~autoParallelFileOperation();


    // Member Functions

        // (reg)IOobject functinality

            //- Search for an object. checkGlobal : also check undecomposed case
            //  This differs from uncollated in searching for parent
            //  also for files in time directories
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&,
                const word& typeName
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "autoParallelFileOperationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
