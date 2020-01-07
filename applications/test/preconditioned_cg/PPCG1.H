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

Class
    Foam::PPCG1

Description
    Preconditioned pipelined conjugate gradient solver for symmetric
    lduMatrices using a run-time selectable preconditioner.

See
    "Hiding global synchronization latency in the preconditioned Conjugate
     Gradient algorithm", P. Ghysels, W. Vanroose

SourceFiles
    PPCG1.C

\*---------------------------------------------------------------------------*/

#ifndef PPCG1_H
#define PPCG1_H

#include "lduMatrix.H"
#include <mpi.h>

#if defined(WM_SP)
    #define MPI_SCALAR MPI_FLOAT
#elif defined(WM_DP)
    #define MPI_SCALAR MPI_DOUBLE
#elif defined(WM_LP)
    #define MPI_SCALAR MPI_LONG_DOUBLE
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PPCG1 Declaration
\*---------------------------------------------------------------------------*/

class PPCG1
:
    public lduMatrix::solver
{
    // Private Member Functions

        void calcDirections
        (
            FixedList<scalar, 3>& globalSum,
            const scalarField& r,
            const scalarField& u,
            const scalarField& w,
            MPI_Request& outstandingRequest
        ) const;

        //- Disallow default bitwise copy construct
        PPCG1(const PPCG1&);

        //- Disallow default bitwise assignment
        void operator=(const PPCG1&);


public:

    //- Runtime type information
    TypeName("PPCG1");


    // Constructors

        //- Construct from matrix components and solver controls
        PPCG1
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~PPCG1()
    {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
