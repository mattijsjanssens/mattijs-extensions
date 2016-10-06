/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::pointPatch

Description
    Basic pointPatch represents a set of points from the mesh.

SourceFiles
    pointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatch_H
#define pointPatch_H

#include "labelList.H"
#include "vectorField.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class pointBoundaryMesh;
class pointConstraint;
class PstreamBuffers;
class pointPatch;

Ostream& operator<<(Ostream&, const pointPatch&);

/*---------------------------------------------------------------------------*\
                      Class pointPatch Declaration
\*---------------------------------------------------------------------------*/

class pointPatch
{
    // Private data

        //- Reference to boundary mesh
        const pointBoundaryMesh& boundaryMesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointPatch(const pointPatch&);

        //- Disallow default bitwise assignment
        void operator=(const pointPatch&);


protected:

    // Protected Member Functions

        // The pointPatch geometry initialisation is called by pointBoundaryMesh
        friend class pointBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&)
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&)
        {}

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&)
        {}


public:

    //- Runtime type information
    TypeName("pointPatch");

    //- Debug switch to disallow the use of genericPointPatch
    static int disallowGenericPointPatch;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatch,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const pointBoundaryMesh& bm,
                const word& patchType
            ),
            (name, dict, index, bm, patchType)
        );


    // Constructor

        pointPatch
        (
            const pointBoundaryMesh& bm
        )
        :
            boundaryMesh_(bm)
        {}


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<pointPatch> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const pointBoundaryMesh&
        );


    //- Destructor
    virtual ~pointPatch()
    {}


    // Member Functions

        //- Return name
        virtual const word& name() const = 0;

        //- Return size
        virtual label size() const = 0;

        //- Return the index of this patch in the pointBoundaryMesh
        virtual label index() const = 0;

        //- Return boundaryMesh reference
        const pointBoundaryMesh& boundaryMesh() const
        {
            return boundaryMesh_;
        }

        //- Return true if this patch field is coupled
        virtual bool coupled() const
        {
            return false;
        }

        //- Return mesh points
        virtual const labelUList& meshPoints() const = 0;

        //- Return mesh points
        virtual const vectorField& localPoints() const = 0;

        //- Return  point normals
        virtual const vectorField& pointNormals() const = 0;

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return word::null;
        }

        //- Accumulate the effect of constraint direction of this patch
        virtual void applyConstraint
        (
            const label pointi,
            pointConstraint&
        ) const
        {}

        //- Write the pointPatch data as a dictionary
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const pointPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
