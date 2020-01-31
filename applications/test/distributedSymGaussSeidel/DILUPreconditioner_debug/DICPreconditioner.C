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

\*---------------------------------------------------------------------------*/

#include "DICPreconditioner.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DICPreconditioner_debug, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<DICPreconditioner_debug>
        addDICPreconditioner_debugSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DICPreconditioner_debug::DICPreconditioner_debug
(
    const lduMatrix::solver& sol,
    const dictionary&
)
:
    lduMatrix::preconditioner(sol),
    rD_(sol.matrix().diag())
{
    const lduInterfaceFieldPtrsList& interfaces = sol.interfaces();
    const FieldField<Field, scalar>& interfaceBouCoeffs =
        sol.interfaceBouCoeffs();

    // Pre-calculate 'left' and 'right' processor interfaces
    DynamicList<label> lowerInterfaces(interfaces.size());
    DynamicList<label> upperInterfaces(interfaces.size());
    {
        forAll(interfaces, inti)
        {
            if
            (
                interfaces.set(inti)
             && isA<processorLduInterface>(interfaces[inti].interface())
            )
            {
                const processorLduInterface& procInt =
                refCast<const processorLduInterface>
                (
                    interfaces[inti].interface()
                );

                Pout<< "    interface:" << inti
                    << " type:" << interfaces[inti].interface().type()
                    << " myProcNo:" << procInt.myProcNo()
                    << " neighbProcNo:" << procInt.neighbProcNo()
                    << endl;

                if (procInt.neighbProcNo() > procInt.myProcNo())
                {
                    upperInterfaces.append(inti);
                }
                else
                {
                    lowerInterfaces.append(inti);
                }
            }
        }
        Pout<< "lowerInterfaces:" << lowerInterfaces << endl;
        Pout<< "upperInterfaces:" << upperInterfaces << endl;
    }


    calcReciprocalD
    (
        rD_,
        sol.matrix(),
        interfaces,
        interfaceBouCoeffs,
        lowerInterfaces,
        upperInterfaces,
        direction(0)                   //cmpt
    );


    // Calculate coefficients for only including lower/upper processors
    
    lowerCoeffs_.setSize(interfaceBouCoeffs.size());
    upperCoeffs_.setSize(interfaceBouCoeffs.size());
    for (const label inti : lowerInterfaces)
    {
        lowerCoeffs_.set
        (
            inti,
            new scalarField(-1.0*interfaceBouCoeffs[inti])
        );
        lowerCoeffs_[inti] *=
            scalarField
            (
                rD_,
                interfaces[inti].interface().faceCells()
            );

        upperCoeffs_.set
        (
            inti,
            new scalarField(-1.0*interfaceBouCoeffs[inti])
        );
        upperCoeffs_[inti] *=
            scalarField
            (
                rD_,
                interfaces[inti].interface().faceCells()
            );                        
    }
    for (const label inti : upperInterfaces)
    {
        lowerCoeffs_.set
        (
            inti,
            new scalarField(interfaceBouCoeffs[inti].size(), 0.0)
        );

        upperCoeffs_.set
        (
            inti,
            new scalarField(interfaceBouCoeffs[inti].size(), 0.0)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DICPreconditioner_debug::calcReciprocalD
(
    scalarField& rD,
    const lduMatrix& matrix,
    const lduInterfaceFieldPtrsList& interfaces,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const labelList& lowerInterfaces,
    const labelList& upperInterfaces,
    const direction cmpt
)
{
    scalar* __restrict__ rDPtr = rD.begin();

    const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();

    // Subtract coupled contributions to the DIC diagonal
    {
        scalarField invRd(matrix.diag());
        forAll(invRd, cell)
        {
            invRd[cell] = 1.0/invRd[cell];
        }

        FieldField<Field, scalar> coeffs(interfaceBouCoeffs.size());
        for (const label inti : lowerInterfaces)
        {
            coeffs.set
            (
                inti,
                new scalarField
                (
                    1.0
                   *interfaceBouCoeffs[inti]
                   *interfaceBouCoeffs[inti]
                )
            );
        }
        for (const label inti : upperInterfaces)
        {
            coeffs.set
            (
                inti,
                new scalarField(interfaceBouCoeffs[inti].size(), 0.0)
            );
        }

        const label startRequest = Pstream::nRequests();
        matrix.initMatrixInterfaces
        (
            coeffs,
            interfaces,
            invRd,
            rD,
            cmpt                
        );
        matrix.updateMatrixInterfaces
        (
            coeffs,
            interfaces,
            invRd,
            rD,
            cmpt,
            startRequest              
        );
    }


    // Calculate the DIC diagonal
    const label nFaces = matrix.upper().size();
    for (label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -= upperPtr[face]*upperPtr[face]/rDPtr[lPtr[face]];
    }


    // Calculate the reciprocal of the preconditioned diagonal
    const label nCells = rD.size();

    for (label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }
}


void Foam::DICPreconditioner_debug::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction cmpt
) const
{
    scalar* __restrict__ wAPtr = wA.begin();
    const scalar* __restrict__ rAPtr = rA.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        solver_.matrix().lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        solver_.matrix().lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr =
        solver_.matrix().upper().begin();

    const lduInterfaceFieldPtrsList& interfaces = solver_.interfaces();

    label nCells = wA.size();
    label nFaces = solver_.matrix().upper().size();
    label nFacesM1 = nFaces - 1;

    for (label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }

    // Do contributions from lower processors
    {
        const label startRequest = Pstream::nRequests();
        solver_.matrix().initMatrixInterfaces
        (
            lowerCoeffs_,
            interfaces,
            wA,
            wA,
            cmpt                
        );
        solver_.matrix().updateMatrixInterfaces
        (
            lowerCoeffs_,
            interfaces,
            wA,
            wA,
            cmpt,
            startRequest              
        );
    }

    for (label face=0; face<nFaces; face++)
    {
        wAPtr[uPtr[face]] -= rDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
    }

    // Do contributions from higher processors
    {
        const label startRequest = Pstream::nRequests();
        solver_.matrix().initMatrixInterfaces
        (
            upperCoeffs_,
            interfaces,
            wA,
            wA,
            cmpt                
        );
        solver_.matrix().updateMatrixInterfaces
        (
            upperCoeffs_,
            interfaces,
            wA,
            wA,
            cmpt,
            startRequest              
        );
    }

    for (label face=nFacesM1; face>=0; face--)
    {
        wAPtr[lPtr[face]] -= rDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
    }
}


// ************************************************************************* //
