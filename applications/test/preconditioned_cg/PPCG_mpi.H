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
    Foam::PPCG

Description
    Preconditioned pipelined conjugate gradient solver for symmetric
    lduMatrices using a run-time selectable preconditioner.

    Reference:
    \verbatim
        P. Ghysels, W. Vanroose.
        "Hiding global synchronization latency in the
         preconditioned Conjugate Gradient algorithm"
    \endverbatim

SourceFiles
    PPCG.C

\*---------------------------------------------------------------------------*/

#ifndef PPCG_H
#define PPCG_H

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
                           Class PPCG Declaration
\*---------------------------------------------------------------------------*/

class PPCG
:
    public lduMatrix::solver
{
    // Private Member Functions

        ////- Non-blocking version of sum(a*b), sum(a*c), sum(mag(sumMag))
        void gSumMagProd
        (
            FixedList<scalar, 3>& globalSum,
            const scalarField& a,
            const scalarField& b,
            const scalarField& c,
            const scalarField& sumMag,
            MPI_Request& outstandingRequest
        ) const;
        //- Non-blocking version of sum(a*b), sum(a*c), sum(mag(sumMag))
        //void gSumMagProd
        //(
        //    FixedList<scalar, 3>& globalSum,
        //    const scalarField& a,
        //    const scalarField& b,
        //    const scalarField& c,
        //    const scalarField& sumMag,
        //    label& outstandingRequest,
        //    const label comm
        //) const;

        //- Disallow default bitwise copy construct
        PPCG(const PPCG&);

        //- Disallow default bitwise assignment
        void operator=(const PPCG&);


protected:

        //- CG solver. Operates either in conjugate-gradient mode or
        //  conjugate residual
        solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt,
            const bool cgMode
        ) const;


public:

    //- Runtime type information
    TypeName("PPCG");


    // Constructors

        //- Construct from matrix components and solver controls
        PPCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~PPCG()
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
