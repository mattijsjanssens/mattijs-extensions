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

Class
    Foam::lduInterface

Description
    An abstract base class for implicitly-coupled interfaces
    e.g. processor and cyclic patches.

SourceFiles
    lduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef lduInterface_H
#define lduInterface_H

#include "labelField.H"
#include "typeInfo.H"
#include "Pstream.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                     Class lduInterface Declaration
\*---------------------------------------------------------------------------*/

class lduInterface
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduInterface(const lduInterface&);

        //- Disallow default bitwise assignment
        void operator=(const lduInterface&);


public:

    //- Runtime type information
    TypeName("lduInterface");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            lduInterface,
            dictionary,
            (
                const dictionary& dict,
                const label index
            ),
            (dict, index)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            lduInterface,
            Istream,
            (
                Istream& is
            ),
            (is)
        );


    // Constructors

        //- Construct null
        lduInterface()
        {}

        //- Construct from Istream
        lduInterface(Istream& is);

        //- Construct and return a subsetted clone
        virtual autoPtr<lduInterface> clone
        (
            const label index,
            const labelUList& mapAddressing
        ) const =  0;


    // Selectors

        //- Return a pointer to a new lduInterface created on freestore from
        //  dictionary
        static autoPtr<lduInterface> New
        (
            const dictionary& dict,
            const label index
        );

        //- Return a pointer to a new lduInterface created on freestore from
        //  stream
        static autoPtr<lduInterface> New(const word&, Istream&);


    //- Destructor
    virtual ~lduInterface();


    // Member Functions

        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const = 0;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const = 0;

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
            ) const = 0;


//        // Mapping functions
//
//            //- Map (and resize as needed) from self given a mapping object
//            virtual void  autoMap(const FieldMapper&) = 0;


        // I/O

            //- Write the data
            virtual void write(Ostream& os) const = 0;

            //- Write the data as a dictionary
            //virtual void writeData(Ostream& os) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
