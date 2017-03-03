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

Class
    Foam::fv::WriteMatrixConstraint

Description
    Constrain the field values within a specified region.

Usage
    For example to set the turbulence properties within a porous region:
    \verbatim
    porosityTurbulence
    {
        type            scalarWriteMatrixConstraint;
        active          yes;

        scalarWriteMatrixConstraintCoeffs
        {
            selectionMode   cellZone;
            cellZone        porosity;
            fieldValues
            {
                k           1;
                epsilon     150;
            }
        }
    }
    \endverbatim

See also
    Foam::fvOption

SourceFiles
    WriteMatrixConstraint.C
    writeMatrixConstraints.C

\*---------------------------------------------------------------------------*/

#ifndef WriteMatrixConstraint_H
#define WriteMatrixConstraint_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class WriteMatrixConstraint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WriteMatrixConstraint
:
    public option
{
    // Private member data


public:

    //- Runtime type information
    TypeName("WriteMatrixConstraint");


    // Constructors

        //- Construct from components
        WriteMatrixConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        //- Set value on field
        virtual void constrain(fvMatrix<Type>& eqn, const label fieldi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WriteMatrixConstraint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
