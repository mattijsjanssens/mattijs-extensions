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
    Foam::staticDisplacementMotionSolver

Description
    Mesh motion solver for an fvMesh. Solves the motion on a separate, static
    region to avoid accumulating discretisation errors.

SourceFiles
    staticDisplacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef staticDisplacementMotionSolver_H
#define staticDisplacementMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class staticDisplacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class staticDisplacementMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Region to solve motion on
        mutable autoPtr<fvMesh> solveMesh_;

        //- Actual solver
        mutable autoPtr<motionSolver> solver_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        staticDisplacementMotionSolver
        (
            const staticDisplacementMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const staticDisplacementMotionSolver&);


public:

    //- Runtime type information
    TypeName("staticDisplacement");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        staticDisplacementMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );


    //- Destructor
    ~staticDisplacementMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
