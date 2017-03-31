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

Class
    Foam::fileOperations::autoDecomposingFileOperation

Description
    fileOperation that assumes file operations are local.

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_autoDecomposingFileOperation_H
#define fileOperations_autoDecomposingFileOperation_H

#include "uncollatedFileOperation.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class fvFieldDecomposer;

namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                Class autoDecomposingFileOperation Declaration
\*---------------------------------------------------------------------------*/

class autoDecomposingFileOperation
:
    public uncollatedFileOperation
{
    mutable autoPtr<Time> baseRunTimePtr_;

    mutable autoPtr<fvMesh> baseMeshPtr_;

    mutable autoPtr<labelIOList> cellAddressingPtr_;

    mutable autoPtr<labelIOList> faceAddressingPtr_;

    mutable autoPtr<labelIOList> boundaryAddressingPtr_;

    mutable autoPtr<fvFieldDecomposer> decomposerPtr_;

    // Private Member Functions

        //- Return the base mesh (demand loaded)
        const fvMesh& baseMesh(const Time&) const;

        //- Return the decomposer (demand loaded)
        const fvFieldDecomposer& decomposer(const IOobject&) const;

        //- Read, decompose and write to stream
        template<class GeoField>
        bool decomposeAndWrite
        (
            const IOobject& procIO,
            const IOobject& parentIO,
            const word& type,
            Ostream& os
        ) const;


public:

        //- Runtime type information
        TypeName("autoDecomposing");


    // Constructors

        //- Construct null
        autoDecomposingFileOperation(const bool verbose);


    //- Destructor
    virtual ~autoDecomposingFileOperation();


    // Member Functions

        // (reg)IOobject functinality

            //- Search for an object. checkGlobal : also check undecomposed case
            //  This differs from uncollated in searching for parent
            //  also for files in time directories
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&
            ) const;

            //- Top-level read
            virtual bool read
            (
                regIOobject&,
                const bool masterOnly,
                const IOstream::streamFormat format,
                const word& typeName
            ) const;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state. Default action is to write to
            //  the objectPath using writeData. If !valid the
            //  file does not need to be written (this is used e.g. to
            //  suppress empty local lagrangian data)
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
    #include "autoDecomposingFileOperationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
