/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::DICPreconditioner_debug

Description
    Simplified diagonal-based incomplete Cholesky preconditioner for symmetric
    matrices (symmetric equivalent of DILU).  The reciprocal of the
    preconditioned diagonal is calculated and stored.

SourceFiles
    DICPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef DICPreconditioner_H
#define DICPreconditioner_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DICPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class DICPreconditioner_debug
:
    public lduMatrix::preconditioner
{
    // Private data

        //- The reciprocal preconditioned diagonal
        scalarField rD_;

        //- Boundary coefficients from lower numbered processors
        FieldField<Field, scalar> lowerCoeffs_;

        //- Boundary coefficients from higher numbered processors
        FieldField<Field, scalar> upperCoeffs_;


public:

    //- Runtime type information
    TypeName("DIC_debug");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        DICPreconditioner_debug
        (
            const lduMatrix::solver&,
            const dictionary& solverControlsUnused
        );


    //- Destructor
    virtual ~DICPreconditioner_debug()
    {}


    // Member Functions

        //- Calculate the reciprocal of the preconditioned diagonal
        static void calcReciprocalD
        (
            scalarField& rD,
            const lduMatrix& matrix,
            const lduInterfaceFieldPtrsList& interfaces,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const labelList& lowerInterfaces,
            const labelList& upperInterfaces,
            const direction cmpt
        );

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
