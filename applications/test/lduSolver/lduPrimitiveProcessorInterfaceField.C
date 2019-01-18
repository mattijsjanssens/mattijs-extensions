/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "lduPrimitiveProcessorInterfaceField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(lduPrimitiveProcessorInterfaceField, 0);
}

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::lduPrimitiveProcessorInterfaceField::lduPrimitiveProcessorInterfaceField
(
    const lduPrimitiveProcessorInterface& interface
)
:
    lduInterfaceField(interface),
    interface_(interface),
    sendBuf_(interface.size()),
    receiveBuf_(interface.size()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //


Foam::lduPrimitiveProcessorInterfaceField::
~lduPrimitiveProcessorInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::lduPrimitiveProcessorInterfaceField::ready() const
{
    if
    (
        outstandingSendRequest_ >= 0
     && outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        bool finished = UPstream::finishedRequest(outstandingSendRequest_);
        if (!finished)
        {
            return false;
        }
    }
    outstandingSendRequest_ = -1;

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        bool finished = UPstream::finishedRequest(outstandingRecvRequest_);
        if (!finished)
        {
            return false;
        }
    }
    outstandingRecvRequest_ = -1;

    return true;
}


void Foam::lduPrimitiveProcessorInterfaceField::initInterfaceMatrixUpdate
(
    scalarField&,
    const scalarField& psiInternal,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    const labelList& fc = interface_.faceCells();
    sendBuf_.setSize(fc.size());
    forAll(fc, i)
    {
        sendBuf_[i] = psiInternal[fc[i]];
    }

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorInFunction
                << "On patch " //<< interface_.name()
                << " outstanding request."
                << abort(FatalError);
        }


        receiveBuf_.setSize(sendBuf_.size());
        outstandingRecvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            interface_.neighbProcNo(),
            reinterpret_cast<char*>(receiveBuf_.begin()),
            receiveBuf_.byteSize(),
            interface_.tag(),
            interface_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            interface_.neighbProcNo(),
            reinterpret_cast<const char*>(sendBuf_.begin()),
            sendBuf_.byteSize(),
            interface_.tag(),
            interface_.comm()
        );
    }
    else
    {
        interface_.compressedSend(commsType, sendBuf_);
    }

    const_cast<lduPrimitiveProcessorInterfaceField&>(*this).updatedMatrix() =
        false;
}


void Foam::lduPrimitiveProcessorInterfaceField::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    const labelUList& faceCells = interface_.faceCells();

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(outstandingRecvRequest_);
        }
        // Recv finished so assume sending finished as well.
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;

        // Consume straight from receiveBuf_

        // Transform according to the transformation tensor
        //transformCoupleField(receiveBuf_, cmpt);

        // Multiply the field by coefficients and add into the result
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*receiveBuf_[elemI];
        }
    }
    else
    {
        scalarField pnf
        (
            interface_.compressedReceive<scalar>(commsType, this->size())()
        );

        // Transform according to the transformation tensor
        //transformCoupleField(pnf, cmpt);

        // Multiply the field by coefficients and add into the result
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }

    const_cast<lduPrimitiveProcessorInterfaceField&>(*this).updatedMatrix() =
        true;
}


// ************************************************************************* //
