/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "scalarField.H"
#include "DICSmoother.H"
#include "DICPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DICSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<DICSmoother>
        addDICSmootherSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DICSmoother::DICSmoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    ),
    rD_(matrix_.diag().size())
{
    const scalarField& diag = matrix_.diag();
    forAll(diag, i)
    {
        rD_[i] = diag[i];
    }

    DICPreconditioner::calcReciprocalD(rD_, matrix_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DICSmoother::smooth
(
    solveScalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    const solveScalar* const __restrict__ rDPtr = rD_.begin();
    const scalar* const __restrict__ upperPtr = matrix_.upper().begin();
    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();

    // Temporary storage for the residual
    solveScalarField rA(rD_.size());
    solveScalar* __restrict__ rAPtr = rA.begin();

    //ConstFieldWrapper<solveScalar, scalar> trD(rD_);

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        matrix_.residual
        (
            rA,
            psi,
            source,
            interfaceBouCoeffs_,
            interfaces_,
            cmpt
        );

        forAll(rA, i)
        {
            rA[i] *= rD_[i];
        }

        label nFaces = matrix_.upper().size();
        for (label facei=0; facei<nFaces; facei++)
        {
            label u = uPtr[facei];
            rAPtr[u] -= rDPtr[u]*upperPtr[facei]*rAPtr[lPtr[facei]];
        }

        label nFacesM1 = nFaces - 1;
        for (label facei=nFacesM1; facei>=0; facei--)
        {
            label l = lPtr[facei];
            rAPtr[l] -= rDPtr[l]*upperPtr[facei]*rAPtr[uPtr[facei]];
        }

        psi += rA;
    }
}


// ************************************************************************* //
