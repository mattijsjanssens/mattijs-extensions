/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019,2022 OpenCFD Ltd.
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
    boundaries. Based on 'Parallel Preconditioners' chapter from
    'Iterative Methods for Sparse Linear Systems' by Yousef Saad.

    Leaves out the handling of boundary nodes after internal nodes since
    probably not beneficial (so no overlap of comms and internal calculation)

See Also
    Foam::DICPreconditioner

SourceFiles
    distributedDICPreconditioner2.C

\*---------------------------------------------------------------------------*/

#ifndef distributedDICPreconditioner2_H
#define distributedDICPreconditioner2_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class distributedDICPreconditioner2 Declaration
\*---------------------------------------------------------------------------*/

class distributedDICPreconditioner2
:
    public lduMatrix::preconditioner
{
protected:

    // Protected data

        //- Buffers for sending and receiving data
        mutable FieldField<Field, solveScalar> sendBufs_;
        mutable FieldField<Field, solveScalar> recvBufs_;
        mutable PtrList<FieldField<Field, solveScalar>> colourBufs_;


        //- Work storage (sized only for non-proc coupled interfaces)
        mutable FieldField<Field, solveScalar> coeffs_;

//        //- Interfaces to non-processor coupled interfaces
//        DynamicList<label> global_;

        //- Interfaces to lower coloured processors
        DynamicList<label> lowerNbrs_;

        //- Interfaces to higher coloured processors
        DynamicList<label> higherNbrs_;

        //- Colour/zone per cell
        labelList cellColour_;

        //- Max colour
        label nColours_;

        //- Per colour the interfaces that (might) influence it

            //- Interfaces to non-processor lower coupled interfaces
            List<DynamicList<label>> lowerGlobal_;

            //- Interfaces to non-processor higher coupled interfaces
            List<DynamicList<label>> higherGlobal_;

//            //- All interfaces
//            List<DynamicList<label>> global_;


        //- The reciprocal preconditioned diagonal
        solveScalarField rD_;


    // Private Member Functions


//- Debug: check that coeffs are zero
void check() const;
        //- Helper: zero out coeffs on selected interfaces
        void zeroCoeffs
        (
            const labelList& selectedInterfaces
        ) const;

        //- Determine zoneing/colouring of local mesh
        void meshColour(const lduMesh& lm, labelList& cellColour) const;

        //- Swap matrix values on selected interfaces
        void updateMatrixInterfaces
        (
            const bool add,
            const FieldField<Field, solveScalar>& coupleCoeffs,
            const labelList& selectedInterfaces,
            const solveScalarField& psiif,
            solveScalarField& result,
            const direction cmpt
        ) const;

        //- Start receiving in recvBufs_
        void receive(const labelList& intfs) const;

        //- Start sending sendBufs_
        void send(const labelList& intfs, const solveScalarField& psi) const;

        //- Send (and store in colourBufs_[colouri]) the effect of
        //  doing selectedInterfaces
        void sendGlobal
        (
            const labelList& selectedInterfaces,
            solveScalarField& psi,
            const label colouri
        ) const;

        //- Calculate reciprocal of diagonal
        void calcReciprocalD(solveScalarField& rD) const;


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

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            solveScalarField& wA,
            const solveScalarField& rA,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
