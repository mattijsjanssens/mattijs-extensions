/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "WriteMatrixConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "fvcCellReduce.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::WriteMatrixConstraint<Type>::WriteMatrixConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::WriteMatrixConstraint<Type>::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.lookup("fields") >> fieldNames_;
        applied_ = boolList(fieldNames_.size(), false);
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void Foam::fv::WriteMatrixConstraint<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "WriteMatrixConstraint<"
        << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const FieldType& psi = eqn.psi();

    if (!psi.time().outputTime())
    {
        return;
    }


    const fvMesh& mesh = psi.mesh();

    {
        word name(psi.name() + "_source");

        Info<< type() << " : Writing source to " << name << endl;
        FieldType source
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>("zero", eqn.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchScalarField::typeName
        );

        source.ref().Field<Type>::operator=(eqn.source());
        source.correctBoundaryConditions();
        source.write();
    }

    if (eqn.hasDiag())
    {
        word name(psi.name() + "_diag");

        Info<< type() << " : Writing diag to " << name << endl;
        volScalarField diag
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        diag.ref().Field<scalar>::operator=(eqn.diag());
        diag.correctBoundaryConditions();
        diag.write();
    }

    if (eqn.hasUpper())
    {
        word name(psi.name() + "_upper");

        Info<< type() << " : Writing upper to " << name << endl;
        surfaceScalarField upper
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );
        upper.ref().Field<scalar>::operator=(eqn.upper());
        upper.write();

        volScalarField minUpper
        (
            psi.name() + "_minUpper",
            fvc::cellReduce(upper, minEqOp<scalar>(), GREAT)
        );
        minUpper.write();

        volScalarField maxUpper
        (
            psi.name() + "_maxUpper",
            fvc::cellReduce(upper, maxEqOp<scalar>(), -GREAT)
        );
        maxUpper.write();
    }

    if (eqn.hasLower())
    {
        word name(psi.name() + "_lower");

        Info<< type() << " : Writing lower to " << name << endl;
        surfaceScalarField lower
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );
        lower.ref().Field<scalar>::operator=(eqn.lower());
        lower.write();

        volScalarField minLower
        (
            psi.name() + "_minLower",
            fvc::cellReduce(lower, minEqOp<scalar>(), GREAT)
        );
        minLower.write();

        volScalarField maxLower
        (
            psi.name() + "_maxLower",
            fvc::cellReduce(lower, maxEqOp<scalar>(), -GREAT)
        );
        maxLower.write();
    }
}


// ************************************************************************* //
