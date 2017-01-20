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
    Foam::fileOperation

\*---------------------------------------------------------------------------*/

#ifndef fileOperation_H
#define fileOperation_H

#include "Istream.H"
#include "Ostream.H"
#include "autoPtr.H"
#include "fileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobject;
class regIOobject;

/*---------------------------------------------------------------------------*\
                         Class fileOperation Declaration
\*---------------------------------------------------------------------------*/

class fileOperation
{

public:

    // Public data types

        //- Enumeration for the location of an IOobject
        enum pathType
        {
            NOTFOUND,               // not found
            ABSOLUTE,               // instance is absolute directory
            OBJECT,                 // objectPath exists
            PROCESSORSOBJECT,       // objectPath exists in processors/
            PARENTOBJECT,           // parent of object path
            FINDINSTANCE,           // file found in time directory
            PROCESSORSFINDINSTANCE  // as above but in processors/
        };


        //- Runtime type information
        TypeName("fileOperation");


    //- Static fileOperation
    static autoPtr<fileOperation> fileHandlerPtr_;


    // Constructors

        //- Construct null
        fileOperation();


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fileOperation,
            word,
            (),
            ()
        );


    // Selectors

        //- Select type
        static autoPtr<fileOperation> New(const word& type);


    //- Destructor
    virtual ~fileOperation();


    // Member Functions

        // OSSpecific equivalents

            //- Make directory
            virtual bool mkDir(const fileName&, mode_t=0777) const = 0;

            //- Set the file mode
            virtual bool chMod(const fileName&, const mode_t) const = 0;

            //- Return the file mode
            virtual mode_t mode(const fileName&) const = 0;

            //- Return the file type: DIRECTORY or FILE
            virtual fileName::Type type(const fileName&) const = 0;

            //- Does the name exist (as DIRECTORY or FILE) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkGzip=true
            ) const = 0;

            //- Does the name exist as a DIRECTORY in the file system?
            virtual bool isDir(const fileName&) const = 0;

            //- Does the name exist as a FILE in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkGzip=true
            ) const = 0;

            //- Return size of file
            virtual off_t fileSize(const fileName&) const = 0;

            //- Return time of last file modification
            virtual time_t lastModified(const fileName&) const = 0;

            //- Return time of last file modification
            virtual double highResLastModified(const fileName&) const = 0;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileName::Type=fileName::FILE,
                const bool filtergz=true
            ) const = 0;

            //- Copy, recursively if necessary, the source to the destination
            virtual bool cp(const fileName& src, const fileName& dst) const = 0;

            //- Create a softlink. dst should not exist. Returns true if
            //  successful.
            virtual bool ln(const fileName& src, const fileName& dst) const = 0;

            //- Rename src to dst
            virtual bool mv(const fileName& src, const fileName& dst) const = 0;

            //- Rename to a corresponding backup file
            //  If the backup file already exists, attempt with
            //  "01" .. "99" suffix
            virtual bool mvBak
            (
                const fileName&,
                const std::string& ext = "bak"
            ) const = 0;

            //- Remove a file, returning true if successful otherwise false
            virtual bool rm(const fileName&) const = 0;

            //- Remove a dirctory and its contents
            virtual bool rmDir(const fileName&) const = 0;

//            //- Open a shared library. Return handle to library. Print error
//            //  message if library cannot be loaded (check = true)
//            virtual void* dlOpen
//            (
//                const fileName& lib,
//                const bool check = true
//            ) const = 0;


        // (reg)IOobject functinality

            //- Search for an object. checkGlobal : also check undecomposed case
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&
            ) const = 0;

            //- Generate an Istream for an IOobject. Returns a null stream
            //  if the file is not valid.
            virtual autoPtr<Istream> objectStream(const fileName&) const = 0;

            //- Reads header for regIOobject and returns an Istream
            //  to read the contents.
            virtual autoPtr<Istream> readStream
            (
                regIOobject&,
                const fileName&
            ) const = 0;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state. Default action is to write to
            //  the objectPath using writeData.
            virtual bool writeObject
            (
                const regIOobject&,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED
            ) const;

            //- Generate an Istream that reads a file
            virtual autoPtr<Istream> NewIFstream(const fileName&) const = 0;

            //- Generate an Ostream that writes a file
            virtual autoPtr<Ostream> NewOFstream
            (
                const fileName& pathname,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global declarations: defined in fileOperation.C

//- Get current file handler
const fileOperation& fileHandler();

//- Reset file handler
void fileHandler(autoPtr<fileOperation>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
