/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "distributedSymGaussSeidelSmoother.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributedSymGaussSeidelSmoother, 0);

    lduMatrix::smoother::
    addsymMatrixConstructorToTable<distributedSymGaussSeidelSmoother>
        adddistributedSymGaussSeidelSmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::
    addasymMatrixConstructorToTable<distributedSymGaussSeidelSmoother>
        adddistributedSymGaussSeidelSmootherAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//void Foam::distributedSymGaussSeidelSmoother::receive
//(
//    const boolList& validInterface,
//    const FieldField<Field, scalar>& coupleCoeffs,
//    const lduInterfaceFieldPtrsList& interfaces,
//    const scalarField& psiif,
//    scalarField& result,
//    const direction cmpt,
//    labelList& outstandinRecvRequest
//)
//{
//    forAll(interfaces, interfacei)
//    {
//        if (interfaces.set(interfacei) && validInterface[interfacei])
//        {
//            UPstream::waitRequest(outstandinRecvRequest[interfacei]);


void Foam::distributedSymGaussSeidelSmoother::initMatrixInterfaces
(
    const lduMatrix& matrix,
    const bool add,
    const boolList& validInterface,
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt
)
{
    forAll(interfaces, interfacei)
    {
        if (interfaces.set(interfacei) && validInterface[interfacei])
        {
            interfaces[interfacei].initInterfaceMatrixUpdate
            (
                result,
                add,
                matrix.mesh().lduAddr(),
                interfacei,
                psiif,
                coupleCoeffs[interfacei],
                cmpt,
                Pstream::defaultCommsType
            );
        }
    }
}


void Foam::distributedSymGaussSeidelSmoother::updateMatrixInterfaces
(
    const lduMatrix& matrix,
    const bool add,
    const boolList& validInterface,
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt
)
{
    // Block for all requests and remove storage
    UPstream::waitRequests();

    // Consume
    forAll(interfaces, interfacei)
    {
        if (interfaces.set(interfacei) && validInterface[interfacei])
        {
            // Processor: receive and add to result
            interfaces[interfacei].updateInterfaceMatrix
            (
                result,
                add,
                matrix.mesh().lduAddr(),
                interfacei,
                psiif,
                coupleCoeffs[interfacei],
                cmpt,
                Pstream::defaultCommsType
            );
        }
    }
}


void Foam::distributedSymGaussSeidelSmoother::forward
(
    scalarField& psi,
    const label celli,
    const lduMatrix& matrix_,
    scalarField& bPrime
)
{
    const scalar* const __restrict__ diagPtr = matrix_.diag().begin();
    const scalar* const __restrict__ upperPtr =
        matrix_.upper().begin();
    const scalar* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();

    scalar* __restrict__ psiPtr = psi.begin();
    scalar* __restrict__ bPrimePtr = bPrime.begin();


    // Start and end of this row
    label fStart = ownStartPtr[celli];
    label fEnd = ownStartPtr[celli + 1];

    // Get the accumulated neighbour side
    scalar psii = bPrimePtr[celli];

    // Accumulate the owner product side
    for (label facei=fStart; facei<fEnd; facei++)
    {
        psii -= upperPtr[facei]*psiPtr[uPtr[facei]];
    }

    // Finish current psi
    psii /= diagPtr[celli];

    // Distribute the neighbour side using current psi
    for (label facei=fStart; facei<fEnd; facei++)
    {
        bPrimePtr[uPtr[facei]] -= lowerPtr[facei]*psii;
    }

    psiPtr[celli] = psii;
}


void Foam::distributedSymGaussSeidelSmoother::back
(
    scalarField& psi,
    const label celli,
    const lduMatrix& matrix_,
    scalarField& bPrime
)
{
    const scalar* const __restrict__ diagPtr = matrix_.diag().begin();
    const scalar* const __restrict__ upperPtr =
        matrix_.upper().begin();
    const scalar* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();

    scalar* __restrict__ psiPtr = psi.begin();
    scalar* __restrict__ bPrimePtr = bPrime.begin();


    // Start and end of this row
    label fEnd = ownStartPtr[celli+1];
    label fStart = ownStartPtr[celli];

    // Get the accumulated neighbour side
    scalar psii = bPrimePtr[celli];

    // Accumulate the owner product side
    for (label facei=fStart; facei<fEnd; facei++)
    {
        psii -= upperPtr[facei]*psiPtr[uPtr[facei]];
    }

    // Finish psi for this cell
    psii /= diagPtr[celli];

    // Distribute the neighbour side using psi for this cell
    for (label facei=fStart; facei<fEnd; facei++)
    {
        bPrimePtr[uPtr[facei]] -= lowerPtr[facei]*psii;
    }

    psiPtr[celli] = psii;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributedSymGaussSeidelSmoother::distributedSymGaussSeidelSmoother
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
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distributedSymGaussSeidelSmoother::smooth
(
    const word& fieldName_,
    scalarField& psi,
    const lduMatrix& matrix_,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs_,
    const lduInterfaceFieldPtrsList& interfaces_,
    const direction cmpt,
    const label nSweeps
)
{
    const label nCells = psi.size();

    scalarField bPrime(nCells);


    // Parallel boundary initialisation.  The parallel boundary is treated
    // as an effective jacobi interface in the boundary.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    FieldField<Field, scalar>& mBouCoeffs =
        const_cast<FieldField<Field, scalar>&>
        (
            interfaceBouCoeffs_
        );

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs[patchi].negate();
        }
    }

    // Determine a few sets of interfaces
    const boolList allInterfaces(interfaces_.size(), true);
    boolList isLowerProcInterface(interfaces_.size(), false);
    boolList isHigherProcInterface(interfaces_.size(), false);

    PackedBoolList toHigherProcCell(psi.size(), false);
    PackedBoolList toLowerProcCell(psi.size(), false);

    forAll(interfaces_, inti)
    {
        if (interfaces_.set(inti))
        {
            const lduInterface& intf = interfaces_[inti].interface();

            if (isA<processorLduInterface>(intf))
            {
                const processorLduInterface& procIf =
                    refCast<const processorLduInterface>(intf);
                isHigherProcInterface[inti] =
                    (procIf.neighbProcNo() > procIf.myProcNo());
                isLowerProcInterface[inti] = !isHigherProcInterface[inti];

                if (procIf.neighbProcNo() > procIf.myProcNo())
                {
                    toHigherProcCell.set(intf.faceCells());
                }
                else
                {
                    toLowerProcCell.set(intf.faceCells());
                }
            }
        }
    }


    // Get (negated) boundary contributions
    scalarField bouVals(psi.size(), 0.0);
    {
        initMatrixInterfaces
        (
            mesh,
            false,
            allInterfaces,
            mBouCoeffs,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );

        updateMatrixInterfaces
        (
            allInterfaces,
            mBouCoeffs,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );
    }


    // Get interior cells
    PackedBoolList isBoundaryCell(psi.size(), false);
    forAll(interfaces_, inti)
    {
        if (interfaces_.set(inti))
        {
            isBoundaryCell.set(interfaces_[inti].interface().faceCells());
        }
    }


    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        bPrime = source;

        // Update cells to lower numbered processors
        updateMatrixInterfaces
        (
            isLowerProcInterface,
            mBouCoeffs,
            interfaces_,
            psi,
            bPrime,
            cmpt
        );

        // Forward solve
        for (label celli=0; celli<nCells; celli++)
        {
            forward(psi, celli, matrix_, bPrime);
        }

        // Start sending to higher numbered processors





        // Forward solve on interior nodes
        for (label celli=0; celli<nCells; celli++)
        {
            if (!isBoundaryCell[celli])
            {
                forward(psi, celli, matrix_, bPrime);
            }
        }

        //matrix_.initMatrixInterfaces
        //(
        //    mBouCoeffs,
        //    interfaces_,
        //    psi,
        //    bPrime,
        //    cmpt
        //);
        //
        //matrix_.updateMatrixInterfaces
        //(
        //    mBouCoeffs,
        //    interfaces_,
        //    psi,
        //    bPrime,
        //    cmpt
        //);

        // Wait for lower numbered processor
        for


        for (label celli=nCells-1; celli>=0; celli--)
        {
            back(psi, celli, matrix_, bPrime);
        }
    }

    // Restore interfaceBouCoeffs_
    forAll(mBouCoeffs, patchi)
    {
        if (interfaces_.set(patchi))
        {
            mBouCoeffs[patchi].negate();
        }
    }
}


void Foam::distributedSymGaussSeidelSmoother::smooth
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        fieldName_,
        psi,
        matrix_,
        source,
        interfaceBouCoeffs_,
        interfaces_,
        cmpt,
        nSweeps
    );
}


// ************************************************************************* //
