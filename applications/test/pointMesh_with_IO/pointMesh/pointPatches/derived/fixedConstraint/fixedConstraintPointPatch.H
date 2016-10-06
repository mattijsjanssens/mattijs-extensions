/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::fixedConstraintPointPatch

Description

SourceFiles
    fixedConstraintPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef fixedConstraintPointPatch_H
#define fixedConstraintPointPatch_H

#include "pointPatch.H"
#include "polyPatch.H"
#include "autoPtr.H"
#include "pointConstraint.H"
#include "patchIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fixedConstraintPointPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedConstraintPointPatch
:
    public pointPatch,
    public patchIdentifier
{

protected:

//    const word name_;
//
//    const label index_;

//    labelList meshPoints_;
//
    // Slices of pointMesh
    label start_;

    label size_;

    SubList<label> meshPoints_;


    List<pointConstraint> constraints_;

    mutable autoPtr<pointField> localPointsPtr_;

    mutable autoPtr<vectorField> pointNormalsPtr_;


public:

    //- Runtime type information
    TypeName("fixedConstraint");


    // Constructors

        //- Construct from dictionary
        fixedConstraintPointPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const pointBoundaryMesh& bm,
            const word& patchType
        );


    //- Destructor
    virtual ~fixedConstraintPointPatch()
    {}


    // Member Functions

        //- Return name
        virtual const word& name() const
        {
            return patchIdentifier::name();
        }

        //- Return size
        virtual label size() const
        {
            return size_;
        }

        //- Return the index of this patch in the pointBoundaryMesh
        virtual label index() const
        {
            return patchIdentifier::index();
        }

        //- Return mesh points
        virtual const labelUList& meshPoints() const
        {
            return meshPoints_;
        }

        //- Return constraints
        virtual const List<pointConstraint>& constraints() const
        {
            return constraints_;
        }

        //- Return pointField of points in patch
        virtual const pointField& localPoints() const;

        //- Return point unit normals
        virtual const vectorField& pointNormals() const;

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return typeName;
        }

        //- Accumulate the effect of constraint direction of this patch
        virtual void applyConstraint
        (
            const label pointi,
            pointConstraint&
        ) const;

        virtual void setConstraints(const List<pointConstraint>&);

        //- Write the pointPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
