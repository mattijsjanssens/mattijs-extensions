/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 M. Janssens
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
    sphericalTensorFieldTest

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvCFD.H"
//#include "skewCorrectionVectors.H"
//#include "volFields.H"
//#include "surfaceFields.H"
//#include "pisoControl.H"
#include "globalIndex.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class some_data
{
public:

    labelList labels_;
    scalarList scalars_;

    some_data(const labelList& labels, const scalarList& scalars)
        : labels_(labels), scalars_(scalars)
    {}


    // 1. Visible requests
    // ~~~~~~~~~~~~~~~~~~~

    bool write
    (
        DynamicList<UPstream::Request>& sendReq,
        const label proci,
        const labelRange& range,
        const int tag = UPstream::msgType(),
        const int communicator = UPstream::worldComm
    ) const
    {
        UOPstream::write
        (
            sendReq.emplace_back(),
            proci,
            SubList<label>(labels_, range.size(), range.start()),
            tag,
            communicator
        );

        UOPstream::write
        (
            sendReq.emplace_back(),
            proci,
            SubList<scalar>(scalars_, range.size(), range.start()),
            tag,
            communicator
        );
        return true;
    }
    bool read
    (
        DynamicList<UPstream::Request>& recvReq,
        const label proci,
        const labelRange& range,
        const int tag = UPstream::msgType(),
        const int communicator = UPstream::worldComm
    )
    {
        UIPstream::read
        (
            recvReq.emplace_back(),
            proci,
            SubList<label>(labels_, range.size(), range.start()),
            tag,
            communicator
        );
        UIPstream::read
        (
            recvReq.emplace_back(),
            proci,
            SubList<scalar>(scalars_, range.size(), range.start()),
            tag,
            communicator
        );
        return true;
    }

/*
    // 2. Hidden requests
    // ~~~~~~~~~~~~~~~~~~

    bool write
    (
        const label proci,
        const int tag = UPstream::msgType(),
        const int communicator = UPstream::worldComm
    ) const
    {
        UOPstream::write
        (
            UPstream::commsType::nonBlocking,
            proci,
            labels_,
            tag,
            communicator
        );
        UOPstream::write
        (
            UPstream::commsType::nonBlocking,
            proci,
            scalars_,
            tag,
            communicator
        );
        return true;
    }
    bool read
    (
        const label proci,
        const labelRange& range,
        const int tag = UPstream::msgType(),
        const int communicator = UPstream::worldComm
    )
    {
        UIPstream::read
        (
            UPstream::commsType::nonBlocking,
            proci,
            SubList<label>(labels_, range.size(), range.start()),
            tag,
            communicator
        );
        UIPstream::read
        (
            UPstream::commsType::nonBlocking,
            proci,
            SubList<scalar>(scalars_, range.size(), range.start()),
            tag,
            communicator
        );
        return true;
    }


    // 3. General UOPstream (or OStream&) interface
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Ostream& writeList(Ostream& os) const
    {
        auto* ops = isA<UOPstream>(os);
        if (ops)
        {
            // Bypass virtual dispatch; data never leaves container
            UOPstream::write
            (
                UPstream::commsType::nonBlocking,
                ops->target(),
                labels_,
                ops->tag(),
                ops->comm()
            );
            UOPstream::write
            (
                UPstream::commsType::nonBlocking,
                ops->target(),
                scalars_,
                ops->tag(),
                ops->comm()
            );
        }
        else
        {
            os << labels_ << scalars_;
        }
        return os;
    }
    Istream& readList(const labelRange& range, Istream& is)
    {
        auto* ips = isA<UIPstream>(is);
        if (ips)
        {
            // Bypass virtual dispatch
            UIPstream::read
            (
                UPstream::commsType::nonBlocking,
                ips->target(),
                SubList<label>(labels_, range.size(), range.start()),
                ips->tag(),
                ips->comm()
            );
            UIPstream::read
            (
                UPstream::commsType::nonBlocking,
                ips->target(),
                SubList<scalar>(scalars_, range.size(), range.start()),
                ips->tag(),
                ips->comm()
            );
        }
        else
        {
            is >> labels_ >> scalars_;
        }
        return is;
    }
*/
    // 4. Return contigous data
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    void data
    (
        DynamicList<void*>& bufPtrs,
        DynamicList<std::streamsize>& byteSizes,
        DynamicList<label>& elemSizes
    )
    {
        bufPtrs.append(labels_.data());
        byteSizes.append(labels_.size()*sizeof(label));
        elemSizes.append(sizeof(label));
        bufPtrs.append(scalars_.data());
        byteSizes.append(scalars_.size()*sizeof(scalar));
        elemSizes.append(sizeof(scalar));
    }
};


void gather
(
    some_data& my_data,
    const globalIndex& all_sizes,
    const int tag = UPstream::msgType(),
    const int comm = UPstream::worldComm
)
{
    const label my_size = my_data.labels_.size();

     // Start sending and receiving data
    if (!UPstream::master(comm))
    {
        DynamicList<UPstream::Request> sendReq;
        my_data.write
        (
            sendReq,
            UPstream::masterNo(),
            labelRange(my_size),
            tag,
            comm
        );
    }
    else
    {
        // Make space for all data from all procs
        my_data.labels_.resize(all_sizes.totalSize());
        my_data.scalars_.resize(all_sizes.totalSize());

        List<DynamicList<UPstream::Request>> recvReq(UPstream::nProcs(comm));

        // Start receiving data from all sub-procs
        label offset = all_sizes.localSize();
        for (const int proci : UPstream::subProcs(comm))
        {
            my_data.read
            (
                recvReq[proci],
                proci,
                labelRange(offset, all_sizes.localSize(proci)),
                tag,
                comm
            );
            offset += all_sizes.localSize(proci);
        }

        // Check and consume
        bitSet outstanding(UPstream::nProcs(), true);
        while (outstanding.any())
        {
            for (const int proci : outstanding)
            {
                if (UPstream::finishedRequests(recvReq[proci]))
                {
                    // Finished. Do something with the data...
                    const label offset = all_sizes.offsets()[proci];
                    const label count = all_sizes.localSize(proci);
                    Pout<< "** from:" << proci
                        << ":received:"
                        << SubList<label>(my_data.labels_, count, offset)
                        << endl;        
                    outstanding.unset(proci);
                }
            }
        }
    }
}
void scatter
(
    some_data& my_data,
    const globalIndex& all_sizes,
    const int tag = UPstream::msgType(),
    const int comm = UPstream::worldComm
)
{
    const label my_size = my_data.labels_.size();

    if (!UPstream::master(comm))
    {
        DynamicList<UPstream::Request> recvReq;
        my_data.read
        (
            recvReq,
            UPstream::masterNo(),
            labelRange(my_size),
            tag,
            comm
        );
    }
    else
    {
        DynamicList<UPstream::Request> sendReq;

        // Start sending to all sub-procs
        label offset = all_sizes.localSize();
        for (const int proci : UPstream::subProcs(comm))
        {
            my_data.write
            (
                sendReq,
                proci,
                labelRange(offset, all_sizes.localSize(proci)),
                tag,
                comm
            );
            offset += all_sizes.localSize(proci);
        }

        // Wait. Nothing to do with the data on master proc, so just wait
        // for sends to finish
        UPstream::waitRequests(sendReq);

        // Trucate to local size
        my_data.labels_.resize(all_sizes.localSize());
        my_data.scalars_.resize(all_sizes.localSize());
    }
}


// Use memory chunks
// ~~~~~~~~~~~~~~~~~

//- Posts sends/receives. Appends to requests.
void inplaceGather
(
    some_data& my_data,
    DynamicList<UPstream::Request>& sendReq,        // slaves only
    List<DynamicList<UPstream::Request>>& recvReq,  // master only
    const globalIndex& all_sizes,                   // master only
    const UList<void*>& datas,
    const UList<std::streamsize>& byteSizes,
    const UList<label>& elemSizes,
    const int tag = UPstream::msgType(),
    const int comm = UPstream::worldComm
)
{
    if (!UPstream::master(comm))
    {
        forAll(datas, i)
        {
            // Send complete chunk
            UOPstream::write
            (
                sendReq.emplace_back(),
                UPstream::masterNo(),
                reinterpret_cast<const char*>(datas[i]),
                byteSizes[i],
                tag,
                comm
            );
        }
    }
    else
    {
        recvReq.resize_nocopy(UPstream::nProcs(comm));

        // Start receiving data from all sub-procs
        label offset = all_sizes.localSize();
        for (const int proci : UPstream::subProcs(comm))
        {
            forAll(datas, i)
            {
                const std::streamsize nBytes =
                    elemSizes[i]*all_sizes.localSize(proci);
                char* recvData =
                    reinterpret_cast<char*>(datas[i])+elemSizes[i]*offset;
                UIPstream::read
                (
                    recvReq[proci].emplace_back(),
                    proci,
                    recvData,
                    nBytes,
                    tag,
                    comm
                );
            }

            offset += all_sizes.localSize(proci);
        }
    }
}


void inplaceGather
(
    some_data& my_data,
    const globalIndex& all_sizes,
    const int tag = UPstream::msgType(),
    const int comm = UPstream::worldComm
)
{
    // Make space for all data from all procs
    if (UPstream::master(comm))
    {
        my_data.labels_.resize(all_sizes.totalSize());
        my_data.scalars_.resize(all_sizes.totalSize());
    }

    // Get my data as contiguous memory chunks
    DynamicList<void*> datas;
    DynamicList<std::streamsize> byteSizes;
    DynamicList<label> elemSizes;
    my_data.data(datas, byteSizes, elemSizes);

    // Start all comms, return requests
    DynamicList<UPstream::Request> sendReq;
    List<DynamicList<UPstream::Request>> recvReq;
    inplaceGather
    (
        my_data,
        sendReq,        // slaves only
        recvReq,        // master only
        all_sizes,      // master only

        datas,
        byteSizes,
        elemSizes,

        tag,
        comm
    );
    UPstream::waitRequests(sendReq);
    for (auto& req : recvReq)
    {
        UPstream::waitRequests(req);
    }
}


//- Single chunk
void inplaceScatter
(
    some_data& my_data,
    List<UPstream::Request>& sendReq,           // master only
    UPstream::Request& recvReq,                 // slave only
    const globalIndex& all_sizes,               // master only
    void* data,
    const std::streamsize byteSize,
    const label elemSize,
    const int tag = UPstream::msgType(),
    const int comm = UPstream::worldComm
)
{
    if (!UPstream::master(comm))
    {
        UIPstream::read
        (
            recvReq,
            UPstream::masterNo(),
            reinterpret_cast<char*>(data),
            byteSize,
            tag,
            comm
        );
    }
    else
    {
        sendReq.resize_nocopy(UPstream::nProcs(comm));

        // Start sending to all sub-procs
        label offset = all_sizes.localSize();
        for (const int proci : UPstream::subProcs(comm))
        {
            const std::streamsize nBytes = elemSize*all_sizes.localSize(proci);
            char* sendData = reinterpret_cast<char*>(data)+elemSize*offset;

            UOPstream::write
            (
                sendReq[proci],
                proci,
                sendData,
                nBytes,
                tag,
                comm
            );

            offset += all_sizes.localSize(proci);
        }
    }
}
//- Multiple same-size chunks
void inplaceScatter
(
    some_data& my_data,
    List<DynamicList<UPstream::Request>>& sendReq,  // master only
    DynamicList<UPstream::Request>& recvReq,        // slave only
    const globalIndex& all_sizes,                   // master only
    const UList<void*>& datas,
    const UList<std::streamsize>& byteSizes,
    const UList<label>& elemSizes,
    const int tag = UPstream::msgType(),
    const int comm = UPstream::worldComm
)
{
    if (!UPstream::master(comm))
    {
        forAll(datas, i)
        {
            UIPstream::read
            (
                recvReq.emplace_back(),
                UPstream::masterNo(),
                reinterpret_cast<char*>(datas[i]),
                byteSizes[i],
                tag,
                comm
            );
        }
    }
    else
    {
        sendReq.resize_nocopy(UPstream::nProcs(comm));

        // Start sending to all sub-procs
        label offset = all_sizes.localSize();
        for (const int proci : UPstream::subProcs(comm))
        {
            for (label i = 0; i < datas.size(); ++i)
            {
                const std::streamsize nBytes =
                    elemSizes[i]*all_sizes.localSize(proci);
                char* sendData =
                    reinterpret_cast<char*>(datas[i])+elemSizes[i]*offset;

                UOPstream::write
                (
                    sendReq[proci].emplace_back(),
                    proci,
                    sendData,
                    nBytes,
                    tag,
                    comm
                );
            }
            offset += all_sizes.localSize(proci);
        }
    }
}
//- Multiple same-size chunks
void inplaceScatter
(
    some_data& my_data,
    const globalIndex& all_sizes,       // assume all the same size
    const int tag = UPstream::msgType(),
    const int comm = UPstream::worldComm
)
{
    // Get my data as contiguous memory chunks
    DynamicList<void*> datas;
    DynamicList<std::streamsize> byteSizes;
    DynamicList<label> elemSizes;
    my_data.data(datas, byteSizes, elemSizes);

    List<DynamicList<UPstream::Request>> sendReq;
    DynamicList<UPstream::Request> recvReq;
    inplaceScatter
    (
        my_data,
        sendReq,        // master only
        recvReq,        // slave only
        all_sizes,      // master only
        datas,
        byteSizes,
        elemSizes,
        tag,
        comm
    );

    // Wait. Do consumption here.
    UPstream::waitRequests(recvReq);
    for (auto& req : sendReq)
    {
        UPstream::waitRequests(req);
    }

    // Trucate to local size
    if (UPstream::master(comm))
    {
        my_data.labels_.resize(all_sizes.localSize());
        my_data.scalars_.resize(all_sizes.localSize());
    }
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // Generate some data on all procs
    some_data my_data
    (
        labelList(3, UPstream::myProcNo()),
        scalarList(3, 1.0*UPstream::myProcNo())
    );
    if (UPstream::myProcNo() % 2)
    {
        my_data.labels_.clear();
        my_data.scalars_.clear();        
    }

    // Gather sizes on master
    const label my_size = my_data.labels_.size();
    const globalIndex all_sizes
    (
        globalIndex::gatherOnly{},
        my_size
    );

    Pout<< "** before:" << flatOutput(my_data.labels_) << endl;

    // // In-place gather data from all procs to master
    // gather(my_data, all_sizes);
    // Pout<< "** after gather:" << flatOutput(my_data.labels_) << endl;

    // // In-place scatter data from master to all procs
    // scatter(my_data, all_sizes);
    // Pout<< "** after scatter:" << flatOutput(my_data.labels_) << endl;


    // Now keep comms outside of the gather/scatter functions, and
    // use the data() interface
    inplaceGather(my_data, all_sizes);
    Info<< "** after inplaceGather:" << flatOutput(my_data.labels_) << endl;

    // Change testdata to make sure the data is actually coming from master
    if (!UPstream::master())
    {
        my_data.labels_ = -1;
        my_data.scalars_ = -1.0;
    }

    Pout<< "** before inplaceScatter:" << flatOutput(my_data.labels_) << endl;
    //inplaceScatter(my_data, all_sizes);
    {
        // Get my data as contiguous memory chunks
        DynamicList<void*> datas;
        DynamicList<std::streamsize> byteSizes;
        DynamicList<label> elemSizes;
        my_data.data(datas, byteSizes, elemSizes);

        // Insert comms
        List<DynamicList<UPstream::Request>> sendReqs(UPstream::nProcs());
        DynamicList<UPstream::Request> recvReqs;

        // Work
        List<UPstream::Request> procToSendReq(UPstream::nProcs());

        forAll(datas, i)
        {
            // Or assume all data elements have same size (in number
            // of elements)?
            const globalIndex all_sizes
            (
                globalIndex::gatherOnly{},
                (
                    UPstream::master()
                  ? my_size
                  : byteSizes[i]/elemSizes[i]
                )
            );

            procToSendReq = UPstream::Request();
            UPstream::Request recvReq;
            inplaceScatter
            (
                my_data,
                procToSendReq,  // master only
                recvReq,        // slave only
                all_sizes,      // master only
                datas[i],
                byteSizes[i],
                elemSizes[i]
            );

            // Append to overall
            forAll(procToSendReq, proci)
            {
                if (procToSendReq[proci].good())
                {
                    sendReqs[proci].append(procToSendReq[proci]);
                }
            }
            if (recvReq.good())
            {
                recvReqs.append(recvReq);
            }
        }
        // Wait. Do consumption here.
        UPstream::waitRequests(recvReqs);
        for (auto& req : sendReqs)
        {
            UPstream::waitRequests(req);
        }

        // Trucate to local size
        if (UPstream::master())
        {
            my_data.labels_.resize(all_sizes.localSize());
            my_data.scalars_.resize(all_sizes.localSize());
        }
        Pout<< "** after inplaceScatter:" << flatOutput(my_data.labels_) << endl;
    }

    return 0;
}


// ************************************************************************* //
