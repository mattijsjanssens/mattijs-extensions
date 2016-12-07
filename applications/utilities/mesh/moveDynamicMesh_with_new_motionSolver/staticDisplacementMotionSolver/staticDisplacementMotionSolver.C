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

\*---------------------------------------------------------------------------*/

#include "staticDisplacementMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(staticDisplacementMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        staticDisplacementMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::staticDisplacementMotionSolver::staticDisplacementMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    solveMesh_
    (
        new fvMesh
        (
            Foam::IOobject
            (
                typeName,           //Foam::fvMesh::defaultRegion,
                mesh.time().timeName(),
                mesh.time(),
                Foam::IOobject::MUST_READ
            )
        )
    ),
    solver_
    (
        motionSolver::New
        (
            solveMesh_(),
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",
                    solveMesh_().time().constant(),
                    solveMesh_(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    true        //false
                )
            )
        )
    )
{
    if (!isA<displacementMotionSolver>(solver_()))
    {
        FatalErrorInFunction
            << "Can only use a displacementMotionSolver type for now"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::staticDisplacementMotionSolver::~staticDisplacementMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::staticDisplacementMotionSolver::curPoints() const
{
    return solver_().curPoints();
}


void Foam::staticDisplacementMotionSolver::solve()
{
    displacementMotionSolver& baseSolver =
        const_cast<displacementMotionSolver&>
        (
            refCast<displacementMotionSolver>(solver_())
        );

    // Update bc
    pointDisplacement().boundaryFieldRef().updateCoeffs();

    // Force onto static mesh. Note: cannot use GeometricField== since checks
    // for same mesh.
    baseSolver.pointDisplacement().ref() == pointDisplacement();
    baseSolver.pointDisplacement().boundaryFieldRef() ==
        pointDisplacement().boundaryField();

    // Solve on static mesh
    solver_().solve();

    // Copy pointDisplacement back
    pointDisplacement().ref() == baseSolver.pointDisplacement();
    pointDisplacement().boundaryFieldRef() ==
        baseSolver.pointDisplacement().boundaryField();
}


void Foam::staticDisplacementMotionSolver::movePoints
(
    const pointField& newPoints
)
{
    //solver_().movePoints(newPoints);
}


void Foam::staticDisplacementMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    solver_().updateMesh(mpm);
}


// ************************************************************************* //
