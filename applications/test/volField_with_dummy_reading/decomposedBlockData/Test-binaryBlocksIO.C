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

Application
    Test-volField

Description
    The idea is that argList allocates a fileServer (through command-line
    arguments or environment vars). This has operations to do
    - file existence checking
    - mkDir, rm etc.
    - open an IFstream
    - open an OFstream
    and everywhere instead of constructing an I/OFstream we do a
        fileServer::NewIFstream(io)
        fileServer::NewOFstream(io)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IOPtrList.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "decomposedBlockData.H"
#include <pthread.h>
#include "FIFOStack.H"
#include "threadedOFstream.H"
#include "OFstreamWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class thread_data
{
public:
    int id_;
    pthread_mutex_t& mutex_;
    pthread_cond_t& condition_;
    FIFOStack<regIOobject*>& objects_;

    thread_data
    (
        const int id,
        pthread_mutex_t& mutex,
        pthread_cond_t& condition,
        FIFOStack<regIOobject*>& objects
    )
    :
        id_(id),
        mutex_(mutex),
        condition_(condition),
        objects_(objects)
    {}
};

void* writeFiles(void *threadarg)
{
    thread_data& my_data = *static_cast<thread_data *>(threadarg);
    DebugVar(my_data.id_);

    FIFOStack<regIOobject*>& objects = my_data.objects_;
    pthread_mutex_t& mutex = my_data.mutex_;
    pthread_cond_t& condition = my_data.condition_;
    const int id = my_data.id_;

    while (true)
    {
        Pout<< id << ":" << "starting consumption" << endl;
        // Consume all objects
        while (true)
        {
            regIOobject* io = nullptr;

            pthread_mutex_lock(&mutex);
            if (objects.size())
            {
                Pout<< id << ":" << "popping" << endl;
                io = objects.pop();
                Pout<< id << ":" << "popped " << io->objectPath() << endl;
            }
            pthread_mutex_unlock(&mutex);

            if (io)
            {
                Pout<< id << ":" << "Consuming " << io->objectPath() << endl;
                sleep(1);
                delete io;
            }
            else
            {
                Pout<< id << ":" << "Empty stack. breaking" << endl;
                break;
            }
        }


        // Wait a bit
        //Pout<< "starting sleep" << endl;
        //sleep(1);
        pthread_mutex_lock(&mutex);
        Pout<< id << ":" << "starting wait" << endl;
        pthread_cond_wait(&condition, &mutex);
        if (!objects.size())
        {
            Pout<< id << ":" << "** signalled to stop?" << endl;
            break;
        }
        Pout<< id << ":" << "finished wait" << endl;
        pthread_mutex_unlock(&mutex);
    }

    Pout<< id << ":" << "exiting thread" << endl;
    pthread_exit(nullptr);
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


{
    OFstreamWriter writeServer(100);

    for (label i = 0; i < 100; i++)
    {
        DebugVar(i);

        threadedOFstream os(writeServer, "junk.txt");
        os << "nlabla" << endl;
    }
}
return 0;

    pthread_t writeThread;
    pthread_mutex_t writeMutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t writeCondition = PTHREAD_COND_INITIALIZER;
    FIFOStack<regIOobject*> objects;

    thread_data data(0, writeMutex, writeCondition, objects);
    int rc = pthread_create(&writeThread, nullptr, writeFiles, &data);

    if (rc)
    {
        FatalErrorInFunction << "problem created thread" << exit(FatalError);
    }

    // Create test field
    for (label i = 0; i < 5; i++)
    {
        DebugVar(i);

        volScalarField* fldPtr = new volScalarField
        (
            IOobject
            (
                "myFld",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );
        volScalarField& fld = *fldPtr;
        fld == dimensionedScalar("zero", dimless, 0.0);

        pthread_mutex_lock(&writeMutex);
        Pout<< "pushed field" << endl;
        objects.push(fldPtr);
        // Signal consumer
        Pout<< "signal consumer" << endl;
        pthread_cond_signal(&writeCondition);
        pthread_mutex_unlock(&writeMutex);
        Pout<< "done pushed field" << endl;

        Pout<< "** starting sleep" << endl;
        sleep(1);
    }


    // Signal thread (by having zero size on stack)
    pthread_mutex_lock(&writeMutex);
    Pout<< "signal consumer to END" << endl;
    pthread_cond_signal(&writeCondition);
    pthread_mutex_unlock(&writeMutex);

    void *status;
    pthread_join(writeThread, &status);
    //pthread_cancel(writeThread);

    pthread_mutex_destroy(&writeMutex);
    pthread_cond_destroy(&writeCondition);
    pthread_exit(nullptr);
    return 0;




Pout<< "std::streamoff:" << sizeof(std::streamoff) << endl;
Pout<< "off_t:" << sizeof(off_t) << endl;
Pout<< "label:" << sizeof(label) << endl;

    IOobject io
    (
        "bufs2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    // Write a field
    {
        volScalarField fld
        (
            IOobject
            (
                "myFld",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );
        fld == dimensionedScalar("zero", dimless, 0.0);

        // Test normal writing
        {
            fld.write();
        }


        // Test collating writing
        OStringStream fldStream;
        fldStream << fld;

        string s(fldStream.str());
        //UList<char> localData(const_cast<char*>(s.data()), s.size());
        //autoPtr<OSstream> osPtr;
        //if (Pstream::master())
        //{
        //    osPtr.reset(new OFstream(io.objectPath(), IOstream::BINARY));
        //    io.writeHeader(osPtr());
        //}
        //List<std::streamoff> start;
        //write(osPtr, localData, start);
        //DebugVar(start);
        UList<char> slice(const_cast<char*>(s.data()), label(s.size()));
        List<char> sList(slice);
        decomposedBlockData localData(io, sList.xfer());
        localData.write();
    }
    {
        io.readOpt() = IOobject::MUST_READ;
        decomposedBlockData localData(io);

        string s(localData.begin(), localData.size());
        IStringStream is(s);
        dictionary dict(is);
        volScalarField fld
        (
            IOobject
            (
                "myFld2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dict
        );
        DebugVar(fld);
    }


//
//
//     // Scan for start and size of entries
//     List<char> localData;
//     {
//         IOList<List<char>> bufs(io);
//
//         autoPtr<Istream> isPtr;
//         if (Pstream::master())
//         {
//             isPtr.reset(new IFstream(io.objectPath()));
//             io.readHeader(isPtr());
//         }
//         read(isPtr, localData);
//
//         DebugVar(localData);
//     }
//
//     // Read a dictionary from the stream
//     {
//         string s(localData.begin(), localData.size());
// DebugVar(s);
//         IStringStream is(s);
//         dictionary dict(is);
//         DebugVar(dict);
//         volScalarField fld
//         (
//             IOobject
//             (
//                 "myFld2",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE,
//                 false
//             ),
//             mesh,
//             dict
//         );
//     }

    Pout<< "**end**" << endl;


    return 0;
}


// ************************************************************************* //
