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
    Foam::fileOperations::mixedCollatedFileOperation

Description
    Version of collatedFileOperation with multiple writers

    - non-parallel writing (e.g. decomposePar):
    Search for existing processorXXX and append to that
    - non-parallel reading (e.g. reconstructPar or foamToVTK -case processor0):
    Search for existing processorXXX and use that
    - parallel writing: use write format
    - parallel reading:
    1) First use write format
    2) Search for processorXXX and see which processor has what. Distribute
       internally if needed

    In general cluster nodes:
    - have no local storage
    - but might have separate paths to NFS directory
    so as long as we multiple cores write (to different files!) we should
    gain speed.

See also
    masterUnmixedCollatedFileOperation

SourceFiles
    mixedCollatedFileOperation.C

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_mixedCollatedFileOperation_H
#define fileOperations_mixedCollatedFileOperation_H

#include "collatedFileOperation.H"
//#include "masterUncollatedFileOperation.H"
//#include "OFstreamCollator.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                    Class mixedCollatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class mixedCollatedFileOperation
:
    public collatedFileOperation    //masterUncollatedFileOperation
{
    // Private data

//        //- Threaded writer
//        mutable OFstreamCollator writer_;

//        //- Number of processors (if not running parallel)
//        label nProcs_;

//        //- Actual name of processors dir
//        word processorsDir_;


   // Private Member Functions

        //- Get the master processor number and size of the communicator group
        static labelPair commsGroup(const label nProcs, const label proci);

        //- Is proc a master
        static bool isMaster(const label nProcs, const label proci);

        //- Get the list of processors part of this set
        static labelList subRanks(const label n);


//        //- Append to processors/ file
//        bool appendObject
//        (
//            const regIOobject& io,
//            const fileName& pathName,
//            IOstream::streamFormat fmt,
//            IOstream::versionNumber ver,
//            IOstream::compressionType cmp
//        ) const;


public:

        //- Runtime type information
        TypeName("mixedCollated");


    // Static data

        //- Max size of thread buffer size. This is the overall size of
        //  all files to be written. Starts blocking if not enough size.
        //  Read as float to enable easy specificiation of large sizes.
        static float maxThreadFileBufferSize;


    // Constructors

        //- Construct null
        mixedCollatedFileOperation(const bool verbose);


    //- Destructor
    virtual ~mixedCollatedFileOperation();


    // Member Functions

        // (reg)IOobject functionality

//            //- Generate disk file name for object. Opposite of filePath.
//            virtual fileName objectPath
//            (
//                const IOobject& io,
//                const word& typeName
//            ) const;

//            //- Writes a regIOobject (so header, contents and divider).
//            //  Returns success state.
//            virtual bool writeObject
//            (
//                const regIOobject&,
//                IOstream::streamFormat format=IOstream::ASCII,
//                IOstream::versionNumber version=IOstream::currentVersion,
//                IOstream::compressionType compression=IOstream::UNCOMPRESSED,
//                const bool valid = true
//            ) const;

        // Other

//            //- Actual name of processors dir. Argument used in non-parallel
//            //  mode to find out processor number
//            virtual word processorsDir(const IOobject&) const;

//            //- Actual name of processors dir. Argument used in non-parallel
//            //  mode to find out processor number
//            virtual word processorsDir(const fileName&) const;

//            //- Set number of processor directories/results. Only used in
//            //  decomposePar
//            virtual void setNProcs(const label nProcs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
