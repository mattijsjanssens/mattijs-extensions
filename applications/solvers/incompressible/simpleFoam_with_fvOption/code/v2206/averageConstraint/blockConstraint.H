/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 M. Janssens
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
    Foam::fv::blockConstraint

Group
    grpFvOptionsConstraints

Description

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    blockConstraint1
    {
        // Mandatory entries (unmodifiable)
        type            blockConstraint;
        selectionMode   cellSet;
        cellSet         c0;
        fields          (T);

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type   | Reqd  | Dflt
      type      | Type name: blockConstraint | word   | yes   | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

See also
  - Foam::Function1

SourceFiles
    blockConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef blockConstraint_H
#define blockConstraint_H

#include "cellSetOption.H"
#include "pointList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class blockConstraint Declaration
\*---------------------------------------------------------------------------*/

class blockConstraint
:
    public fv::cellSetOption
{
    // Private Member Functions

        //void print(fvMatrix<scalar>& eqn, const label celli) const;

        //- Remove connections
        void interpolate
        (
            fvMatrix<scalar>& eqn,
            const label
        );


public:

    //- Runtime type information
    TypeName("blockConstraint");


    // Constructors

        //- Construct from components
        blockConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        blockConstraint(const blockConstraint&) = delete;

        //- No copy assignment
        void operator=(const blockConstraint&) = delete;


    //- Destructor
    virtual ~blockConstraint() = default;


    // Member Functions

        //- Constrain the equation
        virtual void constrain(fvMatrix<scalar>& eqn, const label fieldi);

        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
