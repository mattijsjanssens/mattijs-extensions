/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2022 M. Janssens
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
    Foam::DICPCG

Description
    Preconditioned conjugate gradient solver for symmetric lduMatrices
    using a run-time selectable preconditioner.

SourceFiles
    PCG.C

\*---------------------------------------------------------------------------*/

#ifndef PCG_H
#define PCG_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DICPCG Declaration
\*---------------------------------------------------------------------------*/

class DICPCG
:
    public lduMatrix::solver
{
    // Private Member Functions

        void receive
        (
            const labelList& selectedInterfaces,
            FieldField<Field, scalar>& recvBufs,    // Receive buffer
            label& startOfRequests
        ) const;

        void send
        (
            const labelList& selectedInterfaces,
            const solveScalarField& psiInternal,
            FieldField<Field, scalar>& sendBufs,
            label& startOfRequests
        ) const;


//        void calcReciprocalD
//        (
//            solveScalarField& rD,
//            const lduMatrix& matrix
//        ) const;

        void calcReciprocalD
        (
            solveScalarField& rD,
            const lduMatrix& matrix,
            const labelList& lowerInterfaces,
            const labelList& higherInterfaces,
            FieldField<Field, scalar>& sendBufs,    // Send buffer
            FieldField<Field, scalar>& recvBufs,    // Receive buffer
            label& startOfRequests
        ) const;

        //- Disallow default bitwise copy construct
        DICPCG(const DICPCG&);

        //- Disallow default bitwise assignment
        void operator=(const DICPCG&);


public:

    //- Runtime type information
    TypeName("DICPCG");


    // Constructors

        //- Construct from matrix components and solver controls
        DICPCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~DICPCG()
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
