/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::dummyMotionSolver

Description
    Rigid-body mesh motion solver for fvMesh.

    Applies septernion interpolation of movement as function of distance to the
    object surface.

SourceFiles
    dummyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef dummyMotionSolver_H
#define dummyMotionSolver_H

#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class dummyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class dummyMotionSolver
:
    public motionSolver
{
    // Private Member Functions

        //- No copy construct
        dummyMotionSolver
        (
            const dummyMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const dummyMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("dummyMotionSolver");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        dummyMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        )
        :
            motionSolver(mesh, dict, typeName)
        {}


    //- Destructor
    ~dummyMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const
        {
            return mesh().points();
        }

        //- Solve for motion
        virtual void solve()
        {}

//        //- Write state using given format, version and compression
//        virtual bool writeObject
//        (
//            IOstream::streamFormat fmt,
//            IOstream::versionNumber ver,
//            IOstream::compressionType cmp,
//            const bool valid
//        ) const;
//
//        //- Read dynamicMeshDict dictionary
//        virtual bool read();
//
        //- Update local data for geometry changes
        virtual void movePoints(const pointField&)
        {}

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
