/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::lduPrimitiveProcessorInterfaceField

Description
     Abstract base class for GAMG agglomerated interface fields.

SourceFiles
    lduPrimitiveProcessorInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveProcessorInterfaceField_H
#define lduPrimitiveProcessorInterfaceField_H

#include "lduInterfaceField.H"
#include "lduPrimitiveProcessorInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class lduPrimitiveProcessorInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveProcessorInterfaceField
:
    public lduInterfaceField
{
    // Private data

         //- Local reference cast into the interface
         const lduPrimitiveProcessorInterface& interface_;


        // Sending and receiving

            //- Scalar send buffer
            mutable Field<scalar> sendBuf_;

            //- Scalar receive buffer
            mutable Field<scalar> receiveBuf_;

            //- Outstanding request
            mutable label outstandingSendRequest_;

            //- Outstanding request
            mutable label outstandingRecvRequest_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduPrimitiveProcessorInterfaceField
        (
            const lduPrimitiveProcessorInterfaceField&
        );

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveProcessorInterfaceField&);


public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        lduPrimitiveProcessorInterfaceField
        (
            const lduPrimitiveProcessorInterface& interface
        );


    //- Destructor
    ~lduPrimitiveProcessorInterfaceField();


    // Member Functions

            virtual label size() const
            {
                return interface_.size();
            }

        // Coupled interface matrix update

            //- Is all data available
            virtual bool ready() const;

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
