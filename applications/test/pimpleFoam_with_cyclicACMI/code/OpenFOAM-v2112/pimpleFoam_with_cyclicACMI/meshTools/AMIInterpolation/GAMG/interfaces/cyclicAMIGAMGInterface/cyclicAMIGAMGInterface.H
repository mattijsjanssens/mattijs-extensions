/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::cyclicAMIGAMGInterface

Description
    GAMG agglomerated cyclic AMI interface.

SourceFiles
    cyclicAMIGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIGAMGInterface_H
#define cyclicAMIGAMGInterface_H

#include "GAMGInterface.H"
#include "cyclicAMILduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicAMIGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIGAMGInterface
:
    public GAMGInterface,
    virtual public cyclicAMILduInterface
{
    // Private data

        //- Reference for the cyclicLduInterface from which this is
        //  agglomerated
        const cyclicAMILduInterface& fineCyclicAMIInterface_;

        //- AMI interface
        PtrList<AMIPatchToPatchInterpolation> amiPtrs_;

        //- Overall number of faces in all valid ami
        mutable label neighbSize_;

    // Private Member Functions

        //- No copy construct
        cyclicAMIGAMGInterface(const cyclicAMIGAMGInterface&) = delete;

        //- No copy assignment
        void operator=(const cyclicAMIGAMGInterface&) = delete;


public:

    //- Runtime type information
    TypeName("cyclicAMI");


    // Constructors

        //- Construct from fine level interface,
        //- local and neighbour restrict addressing
        cyclicAMIGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );


    //- Destructor
    virtual ~cyclicAMIGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Cyclic interface functions

            //- Return neighbour patch numbers
            virtual const labelList& neighbPatchIDs() const
            {
                return fineCyclicAMIInterface_.neighbPatchIDs();
            }

            virtual bool owner(const label index) const
            {
                return fineCyclicAMIInterface_.owner(index);
            }

            virtual const cyclicAMIGAMGInterface& neighbPatch
            (
                const label index
            ) const
            {
                return dynamic_cast<const cyclicAMIGAMGInterface&>
                (
                    coarseInterfaces_
                    [
                        fineCyclicAMIInterface_.neighbPatchIDs()[index]
                    ]
                );
            }

            virtual refPtr<AMIPatchToPatchInterpolation> AMI
            (
                const label index
            ) const
            {
                if (!amiPtrs_.set(index))
                {
                    return nullptr;
                }
                else
                {
                    return amiPtrs_[index];
                }
            }

            //- Total neighbour size (of neighbours with valid AMI only)
            virtual label neighbSize() const;

            //- Indices of neighbours with valid AMI
            virtual const labelList& AMIIndices() const
            {
                return fineCyclicAMIInterface_.AMIIndices();
            }

            //- Index on neighbour of AMI
            virtual label neighbIndex(const label index) const
            {
                return fineCyclicAMIInterface_.neighbIndex(index);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return fineCyclicAMIInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return fineCyclicAMIInterface_.reverseT();
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const
            {
                //TBD. How to serialise the AMI such that we can stream
                // cyclicAMIGAMGInterface.
                NotImplemented;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
