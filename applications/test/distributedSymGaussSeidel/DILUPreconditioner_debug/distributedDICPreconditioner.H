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
    Foam::distributedDICPreconditioner

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
    distributedDICPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef distributedDICPreconditioner_H
#define distributedDICPreconditioner_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class distributedDICPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class distributedDICPreconditioner
:
    public lduMatrix::preconditioner
{
protected:

    // Protected data

        //- Buffers for sending and receiving data
        mutable FieldField<Field, solveScalar> sendBufs_;
        mutable FieldField<Field, solveScalar> recvBufs_;
        mutable label startOfRequests_;


        //- Interfaces to lower coloured processors
        DynamicList<label> lowerNbrs_;

        //- Interfaces to higher coloured processors
        DynamicList<label> higherNbrs_;

        //- The reciprocal preconditioned diagonal
        solveScalarField rD_;


    // Private Member Functions

        //- Start receiving in recvBufs_
        void receive(const labelList& intfs) const;

        //- Start sending sendBufs_
        void send(const labelList& intfs, const solveScalarField& psi) const;

        //- Calculate reciprocal of diagonal
        void calcReciprocalD(solveScalarField& rD) const;


public:

    //- Runtime type information
    TypeName("distributedDIC");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        distributedDICPreconditioner
        (
            const lduMatrix::solver&,
            const dictionary& solverControlsUnused
        );


    //- Destructor
    virtual ~distributedDICPreconditioner();


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
