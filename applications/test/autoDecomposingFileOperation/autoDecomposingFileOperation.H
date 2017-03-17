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

    //mutable autoPtr<fvMesh> procMeshPtr_;

    mutable autoPtr<labelIOList> cellAddressingPtr_;

    mutable autoPtr<labelIOList> faceAddressingPtr_;

    mutable autoPtr<labelIOList> boundaryAddressingPtr_;

    mutable autoPtr<fvFieldDecomposer> decomposerPtr_;

    //- Return the base mesh
    const fvMesh& baseMesh(const Time&) const;

    //- Return the current mesh
    const fvMesh& procMesh(const Time&) const;


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
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&
            ) const;

            //- Reads header for regIOobject and returns an ISstream
            //  to read the contents.
            virtual autoPtr<ISstream> readStream
            (
                regIOobject&,
                const fileName&,
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

#endif

// ************************************************************************* //
