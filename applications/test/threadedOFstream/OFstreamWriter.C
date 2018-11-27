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

\*---------------------------------------------------------------------------*/

#include "OFstreamWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstreamWriter, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::OFstreamWriter::writeFile
(
    const fileName& fName,
    const string& data,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool append
)
{
    if (debug)
    {
        Pout<< "OFstreamWriter : Writing " << data.size()
            << " bytes to " << fName << endl;
    }

    mkDir(fName.path());
    OFstream os(fName, fmt, ver, cmp, append);

    if (!os.good())
    {
        return false;
    }
    os.writeQuoted(data, false);

    return os.good();
}


void* Foam::OFstreamWriter::writeAll(void *threadarg)
{
    OFstreamWriter& handler = *static_cast<OFstreamWriter*>(threadarg);

    // Consume stack
    while (true)
    {
        writeData* ptr = nullptr;

        {
            std::lock_guard<std::mutex> guard(handler.mutex_);
            if (handler.objects_.size())
            {
                ptr = handler.objects_.pop();
            }
        }

        if (!ptr)
        {
            break;
        }
        else
        {
            bool ok = writeFile
            (
                ptr->pathName_,
                ptr->data_,
                ptr->format_,
                ptr->version_,
                ptr->compression_,
                ptr->append_
            );
            if (!ok)
            {
                FatalIOErrorInFunction(ptr->pathName_)
                    << "Failed writing " << ptr->pathName_
                    << exit(FatalIOError);
            }

            delete ptr;
        }
        //sleep(1);
    }

    if (debug)
    {
        Pout<< "OFstreamWriter : Exiting write thread " << endl;
    }

    {
        std::lock_guard<std::mutex> guard(handler.mutex_);
        handler.threadRunning_ = false;
    }

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstreamWriter::OFstreamWriter(const off_t maxBufferSize)
:
    maxBufferSize_(maxBufferSize),
    threadRunning_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstreamWriter::~OFstreamWriter()
{
    if (thread_.valid())
    {
        if (debug)
        {
            Pout<< "~OFstreamWriter : Waiting for write thread" << endl;
        }
        thread_().join();
        thread_.clear();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OFstreamWriter::write
(
    const fileName& fName,
    const string& data,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool append
)
{
    if (maxBufferSize_ > 0)
    {
        while (true)
        {
            // Count files to be written
            off_t totalSize = 0;

            {
                std::lock_guard<std::mutex> guard(mutex_);
                forAllConstIter(FIFOStack<writeData*>, objects_, iter)
                {
                    totalSize += iter()->data_.size();
                }
            }

            if
            (
                totalSize == 0
             || (totalSize+off_t(data.size()) < maxBufferSize_)
            )
            {
                break;
            }

            if (debug)
            {
                std::lock_guard<std::mutex> guard(mutex_);
                Pout<< "OFstreamWriter : Waiting for buffer space."
                    << " Currently in use:" << totalSize
                    << " limit:" << maxBufferSize_
                    << " files:" << objects_.size()
                    << endl;
            }

            sleep(5);
        }

        autoPtr<writeData> fileAndDataPtr
        (
            new writeData(fName, data, fmt, ver, cmp, append)
        );

        {
            std::lock_guard<std::mutex> guard(mutex_);

            // Append to thread buffer
            objects_.push(fileAndDataPtr.ptr());

            // Start thread if not running
            if (!threadRunning_)
            {
                if (thread_.valid())
                {
                    if (debug)
                    {
                        Pout<< "OFstreamWriter : Waiting for write thread"
                            << endl;
                    }
                    thread_().join();
                }

                if (debug)
                {
                    Pout<< "OFstreamWriter : Starting write thread"
                        << endl;
                }
                thread_.reset(new std::thread(writeAll, this));
                threadRunning_ = true;
            }
        }

        return true;
    }
    else
    {
        // Immediate writing
        return writeFile(fName, data, fmt, ver, cmp, append);
    }
}


// ************************************************************************* //
