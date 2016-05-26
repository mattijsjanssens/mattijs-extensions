/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Mattijs Janssens
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

#include "displacementFrozenLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementFrozenLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementFrozenLaplacianFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementFrozenLaplacianFvMotionSolver::
displacementFrozenLaplacianFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolverCore(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellDisplacement",
            pointDisplacement_.dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(NULL),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID(coeffDict().lookup("frozenPointsZone"))
      : -1
    )
{
    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "displacementFrozenLaplacianFvMotionSolver:" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }


    if (io.headerOk())
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "displacementFrozenLaplacianFvMotionSolver :"
                << " Read pointVectorField "
                << io.name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }


//     const surfaceScalarField& diff = diffusivity().operator()();
// 
//     fvVectorMatrix eqn
//     (
//         fvm::laplacian
//         (
//             diff,
//             cellDisplacement_,
//             "laplacian(diffusivity,cellDisplacement)"
//         )
//     );

//     laplacianEqn_.reset
//     (
//         new fvVectorMatrix
//         (
//             fvm::laplacian
//             (
//                 diffusivity().operator()(),
//                 cellDisplacement_,
//                 "laplacian(diffusivity,cellDisplacement)"
//             )
//         )
//     );
// 
//     fvVectorMatrix& TEqn = laplacianEqn_();
// 
//     Pout<< "TEqn.diag():" << TEqn.diag() << endl;
//     if (TEqn.hasLower())
//     {
//         Pout<< "TEqn.lower():" << TEqn.lower() << endl;
//     }
//     if (TEqn.hasUpper())
//     {
//         Pout<< "TEqn.upper():" << TEqn.upper() << endl;
//     }
//     Pout<< "TEqn.source():" << TEqn.source() << endl;
//     forAll(TEqn.internalCoeffs(), i)
//     {
//         if (TEqn.internalCoeffs().set(i))
//         {
//             Pout<< "    patch:" << i
//                 << " internal:" << TEqn.internalCoeffs()[i]
//                 << endl;
//         }
//     }
//     forAll(TEqn.boundaryCoeffs(), i)
//     {
//         if (TEqn.boundaryCoeffs().set(i))
//         {
//             Pout<< "    patch:" << i
//                 << " boundary:" << TEqn.boundaryCoeffs()[i]
//                 << endl;
//         }
//     }


//     source_ = TEqn.source();
// 
// DebugVar(TEqn.internalCoeffs().size());
// 
//     const FieldField<Field, vector>& iCoeffs = TEqn.internalCoeffs();
//     internalCoeffs_.setSize(iCoeffs.size());
//     forAll(iCoeffs, i)
//     {
//         if (iCoeffs.set(i))
//         {
//             internalCoeffs_.set(i, iCoeffs[i].clone());
//         }
//     }
// 
// DebugVar(internalCoeffs_.size());
// DebugVar(internalCoeffs_);
// 
// 
//     const FieldField<Field, vector>& bCoeffs = TEqn.boundaryCoeffs();
//     boundaryCoeffs_.setSize(bCoeffs.size());
//     forAll(bCoeffs, i)
//     {
//         if (bCoeffs.set(i))
//         {
//             boundaryCoeffs_.set(i, bCoeffs[i].clone());
//         }
//     }
// 
// DebugVar(boundaryCoeffs_.size());
// DebugVar(boundaryCoeffs_);
// 
// 
// DebugVar("one");
// 
//     faceFluxCorrectionPtr_ = NULL;
//     if (TEqn.faceFluxCorrectionPtr())
//     {
//         faceFluxCorrectionPtr_ =
//             new surfaceVectorField(*TEqn.faceFluxCorrectionPtr());
//     }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementFrozenLaplacianFvMotionSolver::
~displacementFrozenLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionDiffusivity&
Foam::displacementFrozenLaplacianFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_.valid())
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity")
        );
    }
    return diffusivityPtr_();
}


Foam::tmp<Foam::pointField>
Foam::displacementFrozenLaplacianFvMotionSolver::curPoints() const
{
    volPointInterpolation::New(fvMesh_).interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    if (pointLocation_.valid())
    {
        if (debug)
        {
            Info<< "displacementFrozenLaplacianFvMotionSolver : applying "
                << " boundary conditions on " << pointLocation_().name()
                << " to new point location."
                << endl;
        }

        pointLocation_().primitiveFieldRef() =
            points0()
          + pointDisplacement_.primitiveField();

        pointLocation_().correctBoundaryConditions();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                pointLocation_()[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(pointLocation_().primitiveFieldRef());

        return tmp<pointField>(pointLocation_().primitiveField());
    }
    else
    {
        tmp<pointField> tcurPoints
        (
            points0() + pointDisplacement_.primitiveField()
        );
        pointField& curPoints = tcurPoints.ref();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                curPoints[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(curPoints);

        return tcurPoints;
    }
}


void Foam::displacementFrozenLaplacianFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivity().correct();
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

if (!laplacianEqn_.valid())
{
    Pout<< "discretising laplacian" << endl;
    laplacianEqn_.reset
    (
        new fvVectorMatrix
        (
            fvm::laplacian
            (
                diffusivity().operator()(),
                cellDisplacement_,
                "laplacian(diffusivity,cellDisplacement)"
            )
        )
    );
}
fvVectorMatrix eqn
(
    fvm::laplacian
    (
        diffusivity().operator()(),
        cellDisplacement_,
        "laplacian(diffusivity,cellDisplacement)"
    )
);

laplacianEqn_().source() = eqn.source();

Pout<< "solving laplacian" << endl;
Foam::solve(laplacianEqn_());

//    Foam::solve(TEqn);

// 
//     fvVectorMatrix TEqn(cellDisplacement_, cellDisplacement_.dimensions());
// 
//     TEqn.source() = source_;
// 
//     FieldField<Field, vector>& iCoeffs = TEqn.internalCoeffs();
//     iCoeffs.setSize(internalCoeffs_.size());
//     forAll(internalCoeffs_, i)
//     {
//         if (internalCoeffs_.set(i))
//         {
//             iCoeffs[i] = internalCoeffs_[i];
//         }
//     }
// 
//     FieldField<Field, vector>& bCoeffs = TEqn.boundaryCoeffs();
//     bCoeffs.setSize(boundaryCoeffs_.size());
//     forAll(boundaryCoeffs_, i)
//     {
//         if (boundaryCoeffs_.set(i))
//         {
//             bCoeffs[i] = boundaryCoeffs_[i];
//         }
//     }
// 
//     if (faceFluxCorrectionPtr_)
//     {
//         TEqn.faceFluxCorrectionPtr() = *faceFluxCorrectionPtr_;
//     }
//     Foam::solve(laplacianEqn_());
}


void Foam::displacementFrozenLaplacianFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.clear();
}


// ************************************************************************* //
