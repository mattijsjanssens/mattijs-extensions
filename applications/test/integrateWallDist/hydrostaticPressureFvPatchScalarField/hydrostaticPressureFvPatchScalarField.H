/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::hydrostaticPressureFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a hydrostatic pressure condition,
    calculated as:

        \f[
            p_{hyd} = p_{ref} + \rho g (x - x_{ref})
        \f]

    where
    \vartable
        p_{hyd} | Hyrostatic pressure [Pa]
        p_{ref} | Reference pressure [Pa]
        x_{ref} | Reference point in Cartesian co-ordinates
        \rho    | Density
        g       | Acceleration due to gravity [m/s2]
    \endtable

Usage
    \table
        Property     | Description             | Required    | Default value
        rho          | density [kg/m3]         | yes         |
        pRefValue    | reference pressure [Pa] | yes         |
        pRefPoint    | reference pressure location | yes     |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            hydrostaticPressure;
        rho             rho;
        pRefValue       0;
        pRefPoint       (0 0 0);
        value           uniform 0;  // initial value needed since eval needs 'g'
    }
    \endverbatim

Note
    - the pRefPoint is used to find the nearest edge and it adjust
    the pressures such that this edge's value is the reference pressure. There
    is no interpolation if the pRefPoint is not exactly on an edge.
    - it assumes that all of the patch can be reached from a single edge. This
    might not be enough for complex meshes with split outlet patches.

SourceFiles
    hydrostaticPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticPressureFvPatchScalarField_H
#define hydrostaticPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class hydrostaticPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class hydrostaticPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the density field
        word rhoName_;

        //- Reference pressure
        scalar pRefValue_;

        //- Reference pressure location
        vector pRefPoint_;


        // Patch specific information - should be stored on patch

            //- Map from patch edge to coupled edge
            Map<label> patchToCoupledEdge_;

            //- gh on faces
            autoPtr<scalarField> faceCentresPtr_;

            //- gh on edges
            autoPtr<scalarField> edgeCentresPtr_;

            //- Starting edge(s)
            labelList startEdges_;

            //- Edge nearest pRefPoint
            label pRefEdgei_;


    // Private Member Functions

        //- Interpolate onto edges
        tmp<scalarField> linearInterpolate
        (
            const scalarField& faceCentres,
            const scalarField& edgeCentres,
            const scalarField& faceFld
        ) const;

        //- Update topology related properties
        void calcTopology();

        //- Update geometry related properties
        void calcGeometry(const vector& g);


public:

    //- Runtime type information
    TypeName("hydrostaticPressure");


    // Constructors

        //- Construct from patch and internal field
        hydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        hydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  hydrostaticPressureFvPatchScalarField onto a new patch
        hydrostaticPressureFvPatchScalarField
        (
            const hydrostaticPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        hydrostaticPressureFvPatchScalarField
        (
            const hydrostaticPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new hydrostaticPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        hydrostaticPressureFvPatchScalarField
        (
            const hydrostaticPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new hydrostaticPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the reference pressure
            scalar pRefValue() const
            {
                return pRefValue_;
            }

            //- Return reference to the reference pressure to allow adjustment
            scalar& pRefValue()
            {
                return pRefValue_;
            }

            //- Return the pressure reference location
            const vector& pRefPoint() const
            {
                return pRefPoint_;
            }

            //- Return reference to the pressure reference location
            //  to allow adjustment
            vector& pRefPoint()
            {
                return pRefPoint_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
