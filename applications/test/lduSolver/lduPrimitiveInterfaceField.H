/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::lduPrimitiveInterfaceField

Description
    Concrete implementation of lduInterfaceField for procesor transfer.

SourceFiles
    lduPrimitiveInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveInterfaceField_H
#define lduPrimitiveInterfaceField_H

#include "lduPrimitiveInterfaceField.H"
#include "primitiveFieldsFwd.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class lduMatrix;

/*---------------------------------------------------------------------------*\
                     Class lduPrimitiveInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lduPrimitiveInterfaceField
:
    public lduInterfaceField
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduPrimitiveInterfaceField(const lduPrimitiveInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveInterfaceField&);


public:

    //- Runtime type information
    TypeName("lduPrimitiveInterfaceField");


    // Constructors

        //- Construct given coupled patch
        lduPrimitiveInterfaceField(const lduInterface& patch)
        :
            lduInterfaceField(patch)
        {}


    //- Destructor
    virtual ~lduPrimitiveInterfaceField()
    {}


    // Member Functions

        // Coupled interface functionality

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                scalarField&,
                const scalarField&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField&,
                const scalarField&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
