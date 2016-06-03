/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "fixedValuePointPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "pointPatchFieldMapper.H"
#include "pointFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = cf39816173bc28982edc0a1a61d510e46926e5c7
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void rampedFixedValue_cf39816173bc28982edc0a1a61d510e46926e5c7(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchScalarField,
    rampedFixedValueFixedValuePointPatchScalarField
);


const char* const rampedFixedValueFixedValuePointPatchScalarField::SHA1sum =
    "cf39816173bc28982edc0a1a61d510e46926e5c7";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rampedFixedValueFixedValuePointPatchScalarField::
rampedFixedValueFixedValuePointPatchScalarField
(
    const pointPatch& p,
    const DimensionedField<scalar, pointMesh>& iF
)
:
    fixedValuePointPatchField<scalar>(p, iF)
{
    if (#line 30 "/home/mattijs/OpenFOAM/mattijs-extensions/applications/test/codedFunctionObject/cavity/0/pointField.boundaryField.movingWall"
true)
    {
        Info<<"construct rampedFixedValue sha1: cf39816173bc28982edc0a1a61d510e46926e5c7"
            " from patch/DimensionedField\n";
    }
}


rampedFixedValueFixedValuePointPatchScalarField::
rampedFixedValueFixedValuePointPatchScalarField
(
    const rampedFixedValueFixedValuePointPatchScalarField& ptf,
    const pointPatch& p,
    const DimensionedField<scalar, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<scalar>(ptf, p, iF, mapper)
{
    if (#line 30 "/home/mattijs/OpenFOAM/mattijs-extensions/applications/test/codedFunctionObject/cavity/0/pointField.boundaryField.movingWall"
true)
    {
        Info<<"construct rampedFixedValue sha1: cf39816173bc28982edc0a1a61d510e46926e5c7"
            " from patch/DimensionedField/mapper\n";
    }
}


rampedFixedValueFixedValuePointPatchScalarField::
rampedFixedValueFixedValuePointPatchScalarField
(
    const pointPatch& p,
    const DimensionedField<scalar, pointMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValuePointPatchField<scalar>(p, iF, dict, valueRequired)
{
    if (#line 30 "/home/mattijs/OpenFOAM/mattijs-extensions/applications/test/codedFunctionObject/cavity/0/pointField.boundaryField.movingWall"
true)
    {
        Info<<"construct rampedFixedValue sha1: cf39816173bc28982edc0a1a61d510e46926e5c7"
            " from patch/dictionary\n";
    }
}


rampedFixedValueFixedValuePointPatchScalarField::
rampedFixedValueFixedValuePointPatchScalarField
(
    const rampedFixedValueFixedValuePointPatchScalarField& ptf
)
:
    fixedValuePointPatchField<scalar>(ptf)
{
    if (#line 30 "/home/mattijs/OpenFOAM/mattijs-extensions/applications/test/codedFunctionObject/cavity/0/pointField.boundaryField.movingWall"
true)
    {
        Info<<"construct rampedFixedValue sha1: cf39816173bc28982edc0a1a61d510e46926e5c7"
            " as copy\n";
    }
}


rampedFixedValueFixedValuePointPatchScalarField::
rampedFixedValueFixedValuePointPatchScalarField
(
    const rampedFixedValueFixedValuePointPatchScalarField& ptf,
    const DimensionedField<scalar, pointMesh>& iF
)
:
    fixedValuePointPatchField<scalar>(ptf, iF)
{
    if (#line 30 "/home/mattijs/OpenFOAM/mattijs-extensions/applications/test/codedFunctionObject/cavity/0/pointField.boundaryField.movingWall"
true)
    {
        Info<<"construct rampedFixedValue sha1: cf39816173bc28982edc0a1a61d510e46926e5c7 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rampedFixedValueFixedValuePointPatchScalarField::
~rampedFixedValueFixedValuePointPatchScalarField()
{
    if (#line 30 "/home/mattijs/OpenFOAM/mattijs-extensions/applications/test/codedFunctionObject/cavity/0/pointField.boundaryField.movingWall"
true)
    {
        Info<<"destroy rampedFixedValue sha1: cf39816173bc28982edc0a1a61d510e46926e5c7\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rampedFixedValueFixedValuePointPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (#line 30 "/home/mattijs/OpenFOAM/mattijs-extensions/applications/test/codedFunctionObject/cavity/0/pointField.boundaryField.movingWall"
true)
    {
        Info<<"updateCoeffs rampedFixedValue sha1: cf39816173bc28982edc0a1a61d510e46926e5c7\n";
    }

//{{{ begin code
    #line 33 "/home/mattijs/OpenFOAM/mattijs-extensions/applications/test/codedFunctionObject/cavity/0/pointField.boundaryField.movingWall"
operator==
            (
                min(10, 0.1*this->db().time().value())
            );
//}}} end code

    this->fixedValuePointPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

