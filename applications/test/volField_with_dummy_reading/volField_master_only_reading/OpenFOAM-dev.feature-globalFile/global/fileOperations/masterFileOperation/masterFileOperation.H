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
    Foam::fileOperations::masterFileOperation

Description
    fileOperations that performs all file operations on the master processor.
    Requires the calls to be parallel synchronised!

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_masterFileOperation_H
#define fileOperations_masterFileOperation_H

#include "fileOperation.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                     Class masterFileOperation Declaration
\*---------------------------------------------------------------------------*/

class masterFileOperation
:
    public fileOperation
{
protected:

    // Protected classes

        class mkDirOp
        {
            mode_t mode_;
        public:
            mkDirOp(const mode_t mode)
            :
                mode_(mode)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::mkDir(fName, mode_);
            }
        };

        class chModOp
        {
            mode_t mode_;
        public:
            chModOp(const mode_t mode)
            :
                mode_(mode)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::chMod(fName, mode_);
            }
        };

        class modeOp
        {
        public:
            mode_t operator()(const fileName& fName) const
            {
                return Foam::mode(fName);
            }
        };

        class typeOp
        {
        public:
            label operator()(const fileName& fName) const
            {
                return Foam::type(fName);
            }
        };

        class existsOp
        {
            bool checkGzip_;
        public:
            existsOp(const bool checkGzip)
            :
                checkGzip_(checkGzip)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::exists(fName, checkGzip_);
            }
        };

        class isDirOp
        {
        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::isDir(fName);
            }
        };

        class isFileOp
        {
        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::isFile(fName);
            }
        };

        class fileSizeOp
        {
        public:
            off_t operator()(const fileName& fName) const
            {
                return Foam::fileSize(fName);
            }
        };

        class lastModifiedOp
        {
        public:
            time_t operator()(const fileName& fName) const
            {
                return Foam::lastModified(fName);
            }
        };

        class lastModifiedHROp
        {
        public:
            double operator()(const fileName& fName) const
            {
                return Foam::highResLastModified(fName);
            }
        };

        class mvBakOp
        {
            std::string ext_;
        public:
            mvBakOp(const std::string& ext)
            :
                ext_(ext)
            {}

            bool operator()(const fileName& fName) const
            {
                return Foam::mvBak(fName, ext_);
            }
        };

        class rmOp
        {
        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::rm(fName);
            }
        };

        class rmDirOp
        {
        public:
            bool operator()(const fileName& fName) const
            {
                return Foam::rmDir(fName);
            }
        };

        class cpOp
        {
        public:
            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::cp(src, dest);
            }
        };

        class lnOp
        {
        public:
            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::ln(src, dest);
            }
        };

        class mvOp
        {
        public:
            bool operator()(const fileName& src, const fileName& dest) const
            {
                return Foam::mv(src, dest);
            }
        };

        class readDirOp
        {
            const fileName::Type type_;
            const bool filtergz_;
        public:
            readDirOp(const fileName::Type type, const bool filtergz)
            :
                type_(type),
                filtergz_(filtergz)
            {}

            fileNameList operator()(const fileName& fName) const
            {
                return Foam::readDir(fName, type_, filtergz_);
            }
        };


    // Private Member Functions

        template<class Type, class fileOp>
        Type masterOp(const fileName&, const fileOp& fop) const;

        template<class Type, class fileOp>
        Type masterOp
        (
            const fileName&,
            const fileName&,
            const fileOp& fop
        ) const;

        //- Search for filePath; return info on how was found
        static fileName filePath(const bool, const IOobject&, pathType&, word&);

        //- Construct filePath
        static fileName objectPath
        (
            const IOobject&,
            const pathType&,
            const word&
        );


public:

        //- Runtime type information
        TypeName("masterFileOperation");


    // Constructors

        //- Construct null
        masterFileOperation();


    //- Destructor
    virtual ~masterFileOperation();


    // Member Functions

        // OSSpecific equivalents

            //- Make directory
            virtual bool mkDir(const fileName&, mode_t=0777) const;

            //- Set the file mode
            virtual bool chMod(const fileName&, const mode_t) const;

            //- Return the file mode
            virtual mode_t mode(const fileName&) const;

            //- Return the file type: DIRECTORY or FILE
            virtual fileName::Type type(const fileName&) const;

            //- Does the name exist (as DIRECTORY or FILE) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkGzip=true
            ) const;

            //- Does the name exist as a DIRECTORY in the file system?
            virtual bool isDir(const fileName&) const;

            //- Does the name exist as a FILE in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkGzip=true
            ) const;

            //- Return size of file
            virtual off_t fileSize(const fileName&) const;

            //- Return time of last file modification
            virtual time_t lastModified(const fileName&) const;

            //- Return time of last file modification
            virtual double highResLastModified(const fileName&) const;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileName::Type=fileName::FILE,
                const bool filtergz=true
            ) const;

            //- Copy, recursively if necessary, the source to the destination
            virtual bool cp(const fileName& src, const fileName& dst) const;

            //- Create a softlink. dst should not exist. Returns true if
            //  successful.
            virtual bool ln(const fileName& src, const fileName& dst) const;

            //- Rename src to dst
            virtual bool mv(const fileName& src, const fileName& dst) const;

            //- Rename to a corresponding backup file
            //  If the backup file already exists, attempt with
            //  "01" .. "99" suffix
            virtual bool mvBak
            (
                const fileName&,
                const std::string& ext = "bak"
            ) const;

            //- Remove a file, returning true if successful otherwise false
            virtual bool rm(const fileName&) const;

            //- Remove a dirctory and its contents
            virtual bool rmDir(const fileName&) const;

//            //- Open a shared library. Return handle to library. Print error
//            //  message if library cannot be loaded (check = true)
//            virtual void* dlOpen
//            (
//                const fileName& lib,
//                const bool check = true
//            ) const;


        // (reg)IOobject functinality

            //- Search for an object. checkGlobal : also check undecomposed case
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&
            ) const;

            //- Generate an Istream for an IOobject. Returns a null stream
            //  if the file is not valid.
            virtual autoPtr<Istream> objectStream(const fileName&) const;

            //- Reads header for regIOobject and returns an Istream
            //  to read the contents.
            virtual autoPtr<Istream> readStream
            (
                regIOobject&,
                const fileName&
            ) const;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state.
            virtual bool writeObject
            (
                const regIOobject&,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED
            ) const;

            //- Generate an Istream that reads a file
            virtual autoPtr<Istream> NewIFstream(const fileName&) const;

            //- Generate an Ostream that writes a file
            virtual autoPtr<Ostream> NewOFstream
            (
                const fileName& pathname,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED
            ) const;

        // Other

            //- Same file?
            static bool uniformFile(const fileNameList&);

            //- root+casename with any 'processorXXX' replaced by 'processsors'
            static fileName processorsCasePath(const IOobject&);

            //- Like io.path with provided instance and any 'processorXXX'
            //  replaced by 'processsors'
            static fileName processorsPath(const IOobject&, const word&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "masterFileOperationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
