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
    Foam::PPCR

Description
    Preconditioned pipelined conjugate residuals solver for symmetric
    lduMatrices using a run-time selectable preconditioner.

See
    "Hiding global synchronization latency in the preconditioned Conjugate
     Gradient algorithm", P. Ghysels, W. Vanroose

SourceFiles
    PPCR.C

\*---------------------------------------------------------------------------*/

#ifndef PPCR_H
#define PPCR_H

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
                           Class PPCR Declaration
\*---------------------------------------------------------------------------*/

class PPCR
:
    public lduMatrix::solver
{
    // Private Member Functions

        void calcDirections
        (
            FixedList<scalar, 3>& globalSum,
            const scalarField& m,
            const scalarField& r,
            const scalarField& u,
            const scalarField& w,
            MPI_Request& outstandingRequest
        ) const;

        //- Disallow default bitwise copy construct
        PPCR(const PPCR&);

        //- Disallow default bitwise assignment
        void operator=(const PPCR&);


public:

    //- Runtime type information
    TypeName("PPCR");


    // Constructors

        //- Construct from matrix components and solver controls
        PPCR
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~PPCR()
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
