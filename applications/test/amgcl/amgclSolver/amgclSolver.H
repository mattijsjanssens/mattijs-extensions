/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::amgclSolver

Description
    Interface to amgcl.

SourceFiles
    amgclSolver.C

\*---------------------------------------------------------------------------*/

#ifndef amgclSolver_H
#define amgclSolver_H

#include "lduMatrix.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                         Class amgclSolver Declaration
\*---------------------------------------------------------------------------*/

class amgclSolver
:
    public lduMatrix::solver
{
    // Private data

        //- Global numbering for cells
        globalIndex globalNumbering_;

        //- Offsets into crs storage. ptr[celli] is start of data for celli
        std::vector<int> ptr_;


   // Private Member Functions

        //- Convert matrix into CRS format
        void crs
        (
            const scalarField& source,
            std::vector<double>& val,
            std::vector<int>& col,
            std::vector<double>& rhs
        ) const;

public:

    //- Runtime type information
    TypeName("amgcl");


    // Constructors

        //- Construct from matrix components and solver controls
        amgclSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~amgclSolver()
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
