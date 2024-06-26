/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::lduPrimitiveProcessorInterface

Description
    Concrete implementation of processor interface. Used to temporarily store
    settings.

SourceFiles
    lduPrimitiveProcessorInterface.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveProcessorInterface_H
#define lduPrimitiveProcessorInterface_H

//#include "GAMGInterface.H"
#include "lduInterface.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class lduPrimitiveProcessorInterface Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveProcessorInterface
:
    public lduInterface,
    public processorLduInterface
{
    // Private data

        labelList faceCells_;

        //- Communicator to use for parallel communication
        const label comm_;

        //- My processor rank in communicator
        label myProcNo_;

        //- Neighbouring processor rank in communicator
        label neighbProcNo_;

        //- Transformation tensor
        tensorField forwardT_;

        //- Message tag used for sending
        int tag_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduPrimitiveProcessorInterface(const lduPrimitiveProcessorInterface&);

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveProcessorInterface&);


public:

    //- Runtime type information
    TypeName("processor");

    // Constructors

        //- Construct from components
        lduPrimitiveProcessorInterface
        (
            const labelUList& faceCells,
            const label comm,
            const label myProcNo,
            const label neighbProcNo,
            const tensorField& forwardT,
            const int tag
        );

        //- Construct from Istream
        lduPrimitiveProcessorInterface(Istream& is);

        //- Construct from dictionary
        lduPrimitiveProcessorInterface
        (
            const dictionary& dict,
            const label index
        );

        //- Construct given the original patch and a map
        lduPrimitiveProcessorInterface
        (
            const lduPrimitiveProcessorInterface& pp,
            const label index,
            const labelUList& mapAddressing
        );

        //- Construct and return a subsetted clone
        virtual autoPtr<lduInterface> clone
        (
            const label index,
            const labelUList& mapAddressing
        ) const
        {
            return autoPtr<lduInterface>
            (
                new lduPrimitiveProcessorInterface(*this, index, mapAddressing)
            );
        }


    //- Destructor
    virtual ~lduPrimitiveProcessorInterface();


    // Member Functions

        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return faceCells_;
            }

            virtual label size() const
            {
                return faceCells_.size();
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& iF
            ) const
            {
                return tmp<labelField>(new labelField(iF, faceCells_));
            }

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {
                send(commsType, UIndirectList<label>(iF, faceCells_)());
            }

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {
                return receive<label>(commsType, faceCells_.size());
            }


        //- Processor interface functions

            //- Return communicator used for sending
            virtual label comm() const
            {
                return comm_;
            }

            //- Return processor number (rank in communicator)
            virtual int myProcNo() const
            {
                return myProcNo_;
            }

            //- Return neigbour processor number (rank in communicator)
            virtual int neighbProcNo() const
            {
                return neighbProcNo_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return forwardT_;
            }

            //- Return message tag used for sending
            virtual int tag() const
            {
                return tag_;
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
