/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::lduPrimitiveInterface

Description
    Concrete implementation of lduInterface. Used to store one-to-one
    connectivity.

SourceFiles
    lduPrimitiveInterface.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveInterface_H
#define lduPrimitiveInterface_H

#include "lduInterface.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                     Class lduPrimitiveInterface Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveInterface
:
    public lduInterface
{
    // Private data

        labelList faceCells_;

        labelList faceGlobalCells_;

        labelList nbrFaceGlobalCells_;


public:

    //- Runtime type information
    TypeName("lduPrimitiveInterface");


    // Constructors

        //- Construct from components
        lduPrimitiveInterface
        (
            const globalIndex& globalIndex,
            const labelList& fc,
            const labelList& nbrFc
        )
        :
            faceCells_(fc),
            faceGlobalCells_(faceCells_.size()),
            nbrFaceGlobalCells_(nbrFc)
        {
            forAll(faceCells_, i)
            {
                faceGlobalCells_[i] = globalIndex.toGlobal(faceCells_[i]);
            }
        }

        //- Construct from Istream
        lduPrimitiveInterface(Istream& is)
        :
            faceCells_(is),
            faceGlobalCells_(is),
            nbrFaceGlobalCells_(is)
        {}


    //- Destructor
    virtual ~lduPrimitiveInterface()
    {}


    // Member Functions

        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return faceCells_;
            }

            virtual const labelUList& faceGlobalCells() const
            {
                return faceGlobalCells_;
            }

            //- Return faceCell addressing
            virtual const labelUList& nbrFaceGlobalCells() const
            {
                return nbrFaceGlobalCells_;
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& iF
            ) const
            {
                return tmp<labelField>(new labelField(iF, faceCells_));
            }

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {}

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {
                return tmp<labelField>(new labelField(iF, faceCells_));
            }


        // I/O

            virtual Ostream& write(Ostream& os) const
            {
                return os
                    << faceCells_ << token::SPACE
                    << faceGlobalCells_ << token::SPACE
                    << nbrFaceGlobalCells_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
