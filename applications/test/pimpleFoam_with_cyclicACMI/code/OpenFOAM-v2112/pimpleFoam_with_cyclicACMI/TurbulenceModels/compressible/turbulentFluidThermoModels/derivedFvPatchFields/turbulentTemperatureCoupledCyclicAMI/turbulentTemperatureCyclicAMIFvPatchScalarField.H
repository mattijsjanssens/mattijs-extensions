/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::turbulentTemperatureCyclicAMIFvPatchScalarField

Group
    grpCoupledBoundaryConditions

Description

See also
    Foam::cyclicAMIFvPatchField

SourceFiles
    turbulentTemperatureCyclicAMIFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureCyclicAMIFvPatchScalarField_H
#define turbulentTemperatureCyclicAMIFvPatchScalarField_H

#include "cyclicAMIFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "mappedPatchBase.H"
#include "mappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class turbulentTemperatureCyclicAMIFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureCyclicAMIFvPatchScalarField
:
    public cyclicAMIFvPatchField<scalar>,
    public temperatureCoupledBase,
    public mappedPatchBase,
    public mappedPatchFieldBase<scalar>
{
    // Private data

        // From mappedPatchBase

//            //- Region to sample
//            mutable word sampleRegion_;
//
//            //- Patch (if in sampleMode NEARESTPATCH*)
//            mutable word samplePatch_;


        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;
        autoPtr<PatchFunction1<scalar>> thicknessLayer_;

        //- Conductivity of layers
        scalarList kappaLayers_;
        autoPtr<PatchFunction1<scalar>> kappaLayer_;


public:

    //- Runtime type information
    TypeName("compressible:turbulentTemperatureCyclicAMI");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureCyclicAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureCyclicAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCyclicAMIFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureCyclicAMIFvPatchScalarField
        (
            const turbulentTemperatureCyclicAMIFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentTemperatureCyclicAMIFvPatchScalarField
        (
            const turbulentTemperatureCyclicAMIFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        turbulentTemperatureCyclicAMIFvPatchScalarField
        (
            const turbulentTemperatureCyclicAMIFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member functions

        // Access

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return cyclicAMIFvPatchField<scalar>::type();
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


        // Evaluation functions

            //- Local deltaCoeffs (distance to face, not neighbour cell)
        	tmp<scalarField> patchDeltaCoeffs() const;

            tmp<scalarField> patchKappa
            (
                const scalarField& deltaCoeffs,
                const scalarField& Tp
            ) const;

            //- Return patch-normal gradient
            virtual tmp<Field<scalar>> snGrad
            (
                const scalarField& deltaCoeffs
            ) const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<scalar>&,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<scalar>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
