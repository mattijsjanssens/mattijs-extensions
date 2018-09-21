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
    Foam::fileOperations::autoReconstructFileOperation

Description

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_autoReconstructFileOperation_H
#define fileOperations_autoReconstructFileOperation_H

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
                Class autoReconstructFileOperation Declaration
\*---------------------------------------------------------------------------*/

class autoReconstructFileOperation
:
    public uncollatedFileOperation
{
    // Private Member Data

        //- Caching of objects per time. Updated by readObjects
        mutable HashPtrTable<fileNameList, fileName> procObjects_;


    // Private Member Functions

        //- Check for equivalent object inside processor0/
        bool haveProcPath(const IOobject& io, fileName&) const;

        //- Rewrite lagrangian to processor directory
        fileName equivalentLagrangian(const fileName& dir) const;

        //- Helper to write container elements without spaces
        template<class T>
        static Ostream& write(Ostream&, const T&);

public:

        //- Runtime type information
        TypeName("autoReconstruct");


    // Constructors

        //- Construct null
        autoReconstructFileOperation(const bool verbose);


    //- Destructor
    virtual ~autoReconstructFileOperation();


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

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileType = fileType::file,
                const bool filterVariants = true,
                const bool followLink = true
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

            //- Reads header for regIOobject and returns an ISstream
            //  to read the contents.
            virtual autoPtr<ISstream> readStream
            (
                regIOobject&,
                const fileName&,
                const word& typeName,
                const bool valid = true
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

#ifdef NoRepository
    #include "autoReconstructFileOperationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
