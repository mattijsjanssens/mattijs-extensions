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

    bool write
    (
        DynamicList<UPstream::Request>& sendReq,
        const label proci,
        const int tag = UPstream::msgType(),
        const int communicator = UPstream::worldComm
    ) const
    {
        UOPstream::write
        (
            sendReq.emplace_back(),
            proci,
            labels_,
            tag,
            communicator
        );

        UOPstream::write
        (
            sendReq.emplace_back(),
            proci,
            scalars_,
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
};

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
    const label my_size = my_data.labels_.size();

    // Gather sizes
    const globalIndex all_sizes(globalIndex::gatherOnly{}, my_size);

    // Start sending and receiving data
    DynamicList<UPstream::Request> sendReq;
    DynamicList<UPstream::Request> recvReq;
    if (!UPstream::master())
    {
        my_data.write
        (
            sendReq,
            UPstream::masterNo()
        );
    }
    else
    {
        Pout<< "** before:" << flatOutput(my_data.labels_) << endl;

        // Make space for all data from all procs
        my_data.labels_.resize(all_sizes.totalSize());
        my_data.scalars_.resize(all_sizes.totalSize());

        label offset = my_size;
        for (const int proci : UPstream::subProcs())
        {
            my_data.read
            (
                recvReq,
                proci,
                labelRange(offset, all_sizes.sizes()[proci])
            );
            offset += all_sizes.sizes()[proci];
        }

        // Wait for all sends and receives to finish
        // Note: if we want to overlap with comms we need to keep track
        // of the requests per processor!
        UPstream::waitRequests(sendReq);
        UPstream::waitRequests(recvReq);
        Pout<< "** after:" << flatOutput(my_data.labels_) << endl;
    }

    return 0;
}


// ************************************************************************* //
