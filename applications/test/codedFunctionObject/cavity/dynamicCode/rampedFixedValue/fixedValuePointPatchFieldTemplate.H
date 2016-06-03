/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

Description
    Template for use with dynamic code generation of a
    fixedValue pointPatchField.

    - without state

SourceFiles
    fixedValuePointPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValuePointPatchTemplateScalarField_H
#define fixedValuePointPatchTemplateScalarField_H

#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                A templated FixedValuePointPatchField
\*---------------------------------------------------------------------------*/

class rampedFixedValueFixedValuePointPatchScalarField
:
    public fixedValuePointPatchField<scalar>
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("rampedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        rampedFixedValueFixedValuePointPatchScalarField
        (
            const pointPatch&,
            const DimensionedField<scalar, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rampedFixedValueFixedValuePointPatchScalarField
        (
            const pointPatch&,
            const DimensionedField<scalar, pointMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping a copy onto a new patch
        rampedFixedValueFixedValuePointPatchScalarField
        (
            const rampedFixedValueFixedValuePointPatchScalarField&,
            const pointPatch&,
            const DimensionedField<scalar, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        rampedFixedValueFixedValuePointPatchScalarField
        (
            const rampedFixedValueFixedValuePointPatchScalarField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<scalar>> clone() const
        {
            return autoPtr<pointPatchField<scalar>>
            (
                new rampedFixedValueFixedValuePointPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rampedFixedValueFixedValuePointPatchScalarField
        (
            const rampedFixedValueFixedValuePointPatchScalarField&,
            const DimensionedField<scalar, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<scalar>> clone
        (
            const DimensionedField<scalar, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<scalar>>
            (
                new rampedFixedValueFixedValuePointPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~rampedFixedValueFixedValuePointPatchScalarField();


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

