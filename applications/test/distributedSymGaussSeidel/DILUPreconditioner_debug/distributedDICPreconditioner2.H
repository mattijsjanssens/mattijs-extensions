/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::distributedDICPreconditioner2

Group
    grpLduMatrixPreconditioners

Description
    Version of DICpreconditioner that uses preconditioning across processor
    (and coupled) boundaries. Based on 'Parallel Preconditioners' chapter from
    'Iterative Methods for Sparse Linear Systems' by Yousef Saad.

    Leaves out the handling of boundary nodes after internal nodes since
    probably not beneficial (so no overlap of comms and internal calculation)

    By defaults preconditions across coupled boundaries (currently only
    cyclicAMI). This can be disabled with

        solver          PCG;
        preconditioner
        {
            preconditioner  distributedDIC;
            coupled     false;
        }

See Also
    Foam::DICPreconditioner

SourceFiles
    distributedDICPreconditioner2.C

\*---------------------------------------------------------------------------*/

#ifndef distributedDICPreconditioner2_H
#define distributedDICPreconditioner2_H

#include "distributedDILUPreconditioner2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class distributedDICPreconditioner2 Declaration
\*---------------------------------------------------------------------------*/

class distributedDICPreconditioner2
:
    public distributedDILUPreconditioner2
{
protected:

    // Protected data

//        //- Precondition across global coupled bc
//        const bool coupled_;
//
//
//        //- Processor interface buffers and colouring
//
//            //- Buffers for sending and receiving data
//            mutable FieldField<Field, solveScalar> sendBufs_;
//            mutable FieldField<Field, solveScalar> recvBufs_;
//            mutable DynamicList<UPstream::Request> recvRequests_;
//
//            //- Interfaces to lower coloured processors
//            DynamicList<label> lowerNbrs_;
//            mutable DynamicList<UPstream::Request> lowerSendRequests_;
//            mutable DynamicList<UPstream::Request> lowerRecvRequests_;
//
//            //- Interfaces to higher coloured processors
//            DynamicList<label> higherNbrs_;
//            mutable DynamicList<UPstream::Request> higherSendRequests_;
//            mutable DynamicList<UPstream::Request> higherRecvRequests_;
//
//
//        //- Local (cell) colouring from global interfaces
//
//            //- Colour/zone per cell
//            autoPtr<labelList> cellColourPtr_;
//
//            //- Number of colours (in case of multiple disconnected regions
//            //  in single mesh)
//            label nColours_;
//
//
//        //- Global interfaces. Per colour the interfaces that (might)
//        //- influence it
//
//            mutable PtrList<FieldField<Field, solveScalar>> colourBufs_;
//
//            //- Interfaces to non-processor lower coupled interfaces
//            List<DynamicList<label>> lowerGlobalRecv_;
//
//            //- Interfaces to non-processor lower coupled interfaces
//            List<DynamicList<label>> lowerGlobalSend_;
//
//            //- Corresponding destination colour (for lowerGlobal)
//            List<label> lowerColour_;
//
//            //- Interfaces to non-processor higher coupled interfaces
//            List<DynamicList<label>> higherGlobalRecv_;
//
//            //- Interfaces to non-processor higher coupled interfaces
//            List<DynamicList<label>> higherGlobalSend_;
//
//            //- Corresponding destination colour (for higherGlobal)
//            List<label> higherColour_;
//
//
//        //- The reciprocal preconditioned diagonal
//        solveScalarField rD_;
//
//
//    // Private Member Functions
//
//        //- Variant of lduMatrix::updateMatrixInterfaces on selected interfaces
//        void updateMatrixInterfaces
//        (
//            const bool add,
//            const FieldField<Field, solveScalar>& coupleCoeffs,
//            const labelList& selectedInterfaces,
//            const solveScalarField& psiif,
//            solveScalarField& result,
//            const direction cmpt
//        ) const;
//
//        //- Send (and store in colourBufs_[colouri]) the effect of
//        //  doing selectedInterfaces
//        void sendGlobal
//        (
//            const labelList& selectedInterfaces,
//            solveScalarField& psi,
//            const label colouri
//        ) const;
//
//
//        //- Start receiving in recvBufs_
//        void receive
//        (
//            const labelList& selectedInterfaces,
//            DynamicList<UPstream::Request>& requests
//        ) const;
//
//        //- Start sending sendBufs_
//        void send
//        (
//            const labelList& selectedInterfaces,
//            const solveScalarField& psiInternal,
//            DynamicList<UPstream::Request>& requests
//        ) const;
//
//        //- Wait for explicit requests
//        void wait
//        (
//            DynamicList<UPstream::Request>& requests
//        ) const;

        //- Update diagonal for colour
        virtual void forwardInternalDiag
        (
            solveScalarField& rD,
            const label colouri
        ) const;

        //- Update preconditioned variable walking forward on internal faces
        virtual void forwardInternal
        (
            solveScalarField& wA,
            const label colouri
        ) const;

//        //- Update preconditioned variable walking backward on internal faces
//        virtual void backwardInternal
//        (
//            solveScalarField& wA,
//            const label colouri
//        ) const;
//
//        //- Calculate reciprocal of diagonal
//        void calcReciprocalD(solveScalarField& rD) const;


public:

    //- Runtime type information
    TypeName("distributedDIC2");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        distributedDICPreconditioner2
        (
            const lduMatrix::solver&,
            const dictionary& solverControlsUnused
        );


    //- Destructor
    virtual ~distributedDICPreconditioner2() = default;


    // Member Functions

//        //- Return wA the preconditioned form of residual rA
//        virtual void precondition
//        (
//            solveScalarField& wA,
//            const solveScalarField& rA,
//            const direction cmpt=0
//        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
