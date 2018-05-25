/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::fileOperations::autoReconstructingFileOperation

Description

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_autoReconstructingFileOperation_H
#define fileOperations_autoReconstructingFileOperation_H

#include "uncollatedFileOperation.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                Class autoReconstructingFileOperation Declaration
\*---------------------------------------------------------------------------*/

class autoReconstructingFileOperation
:
    public uncollatedFileOperation
{
    // Private Member Data

        //- Caching of objects per time. Updated by readObjects
        mutable HashPtrTable<fileNameList, fileName> procObjects_;


    // Private Member Functions

        //- Check for equivalent object inside processor0/
        bool haveProcPath(const IOobject& io, fileName&) const;

//         //- Reconstruct volume field and stream
//         template<class Type>
//         static Ostream& writeReconstructedFvVolumeField
//         (
//             const fvMesh& mesh,
//             const IOobject& io,
//             Ostream& os
//         );
//
//         //- Reconstruct surface field and stream
//         template<class Type>
//         static Ostream& writeReconstructedFvSurfaceField
//         (
//             const fvMesh& mesh,
//             const IOobject& io,
//             Ostream& os
//         );
//
//
public:

        //- Runtime type information
        TypeName("autoReconstructing");


    // Constructors

        //- Construct null
        autoReconstructingFileOperation(const bool verbose);


    //- Destructor
    virtual ~autoReconstructingFileOperation();


    // Member Functions

        // fileOperation functionality

            //- Reuse file-based filePath
            using fileOperation::filePath;

            //- Search for an object. checkGlobal : also check undecomposed case
            //  This differs from uncollated in searching for parent
            //  also for files in time directories
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&,
                const word& typeName
            ) const;

            //- Search for a directory. checkGlobal : also check undecomposed
            //  case
            virtual fileName dirPath
            (
                const bool checkGlobal,
                const IOobject& io
            ) const;

            //- Get sorted list of times
            virtual instantList findTimes(const fileName&, const word&) const;

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
