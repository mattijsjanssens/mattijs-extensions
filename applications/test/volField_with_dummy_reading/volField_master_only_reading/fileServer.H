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
    Foam::fileServer

\*---------------------------------------------------------------------------*/

#ifndef fileServer_H
#define fileServer_H

#include "Istream.H"
#include "Ostream.H"
#include "autoPtr.H"
#include "IOobject.H"
//#include "edgeList.H"
//#include "fileServerFormatsCore.H"
//#include "runTimeSelectionTables.H"
//#include "memberFunctionSelectionTables.H"
//#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//// Forward declaration of classes
//class Istream;
//class Ostream;
//
//// Forward declaration of friend functions and operators
//class fileServer;
//Istream& operator>>(Istream&, fileServer&);
//Ostream& operator<<(Ostream&, const fileServer&);


/*---------------------------------------------------------------------------*\
                         Class fileServer Declaration
\*---------------------------------------------------------------------------*/

class fileServer
{

public:

    // Public data types

        //- Enumeration for the location of an IOobject
        enum pathType
        {
            NOTFOUND,           // not found
            ABSOLUTE,           // instance is absolute directory
            OBJECT,             // objectPath exists
            PARENTOBJECT,       // parent of object path
            FINDINSTANCE        // file found in time directory
        };


private:
//    // Private data
//
//        //- Vertices of the edges
//        pointField points_;
//
//        //- The edges defining the boundary
//        edgeList edges_;
//
//        //- From point to edges
//        mutable autoPtr<labelListList> pointEdgesPtr_;
//
//
//    // Private Member Functions
//
//        //- Calculate point-edge addressing (inverse of edges)
//        void calcPointEdges() const;
//

public:

        //- Runtime type information
        TypeName("fileServer");


    // Constructors

        //- Construct null
        fileServer();


//    // Declare run-time constructor selection table
//
//        declareRunTimeSelectionTable
//        (
//            autoPtr,
//            fileServer,
//            fileExtension,
//            (),
//            ()
//        );
//
//
//    // Selectors
//
//        //- Select type
//        static autoPtr<fileServer> New
//        (
//            const word& type
//        );


    //- Destructor
    virtual ~fileServer();


    // Member Function Selectors

//        declareMemberFunctionSelectionTable
//        (
//            void,
//            fileServer,
//            write,
//            fileExtension,
//            (
//                const fileName& name,
//                const fileServer& mesh
//            ),
//            (name, mesh)
//        );

    // Member Functions

        //- Generate an Istream that reads a file
        virtual autoPtr<Istream> NewIFstream(IOobject& io) const = 0;

        //- Generate an Ostream that reads a file
        virtual autoPtr<Ostream> NewOFstream
        (
            const fileName& pathname,
            IOstream::streamFormat format=IOstream::ASCII,
            IOstream::versionNumber version=IOstream::currentVersion,
            IOstream::compressionType compression=IOstream::UNCOMPRESSED
        ) const = 0;

        //- Make directory
        virtual bool mkDir(const fileName&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
