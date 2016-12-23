/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "simpleSolver.H"
#include "fvCFD.H"
#include "meshToMesh0.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleSolver::simpleSolver(const fvMesh& mesh)
:
    runTime_(mesh.time()),
    mesh_(mesh),

    simple_(const_cast<fvMesh&>(mesh_)),

    p_
    (
        IOobject
        (
            "p",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    U_
    (
        IOobject
        (
            "U",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    phi_
    (
        IOobject
        (
            "phi",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U_) & mesh_.Sf()
    ),

    pRefCell_(0),
    pRefValue_(0),

    laminarTransport_(U_, phi_),

    turbulence_
    (
        incompressible::turbulenceModel::New(U_, phi_, laminarTransport_)
    ),

    MRF_(mesh_),
    fvOptions_(mesh_)
{
    setRefCell(p_, simple_.dict(), pRefCell_, pRefValue_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleSolver::~simpleSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> Foam::simpleSolver::UEqn(const bool addErr)
{
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi_, U_)
      + MRF_.DDt(U_)
      + turbulence_->divDevReff(U_)
     ==
        fvOptions_(U_)
    );

    UEqn.ref().relax();

    fvOptions_.constrain(UEqn.ref());

    if (addErr && Uerr_.valid())
    {
        UEqn.ref() += Uerr_();
    }

    return UEqn;
}


void Foam::simpleSolver::solve(const int nIter)
{
    MRF_.correctBoundaryVelocity(U_);

    for (int iter=0; iter<nIter; iter++)
    {
        p_.storePrevIter();

        // --- Pressure-velocity SIMPLE corrector
        {
            tmp<fvVectorMatrix> UEqn(this->UEqn());
            ::solve(UEqn() + fvc::grad(p_));
            fvOptions_.correct(U_);

            volScalarField rAU(1.0/UEqn().A());
            volVectorField HbyA("HbyA", U_);
            HbyA = rAU*UEqn().H();

            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::interpolate(HbyA) & mesh_.Sf()
            );
            MRF_.makeRelative(phiHbyA);
            adjustPhi(phiHbyA, U_, p_);

            tmp<volScalarField> rAtU(rAU);

            if (simple_.consistent())
            {
                rAtU = 1.0/(1.0/rAU - UEqn().H1());
                phiHbyA +=
                    fvc::interpolate(rAtU() - rAU)
                   *fvc::snGrad(p_)*mesh_.magSf();
                HbyA -= (rAU - rAtU())*fvc::grad(p_);
            }

            UEqn.clear();

            // Non-orthogonal pressure corrector loop
            while (simple_.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAtU(), p_) == fvc::div(phiHbyA)
                );

                if (pErr_.valid())
                {
                    pEqn -= pErr_();
                }

                pEqn.setReference(pRefCell_, pRefValue_);

                pEqn.solve();

                if (simple_.finalNonOrthogonalIter())
                {
                    phi_ = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p_.relax();

            // Momentum corrector
            U_ = HbyA - rAtU()*fvc::grad(p_);
            U_.correctBoundaryConditions();
            fvOptions_.correct(U_);
        }

        laminarTransport_.correct();
        turbulence_->correct();
    }
}


tmp<volVectorField> Foam::simpleSolver::Ures(const bool addErr)
{
    return (UEqn(addErr) + fvc::grad(p_)) & U_;
}


tmp<volScalarField> Foam::simpleSolver::pRes(const bool addErr)
{
    if (addErr && pErr_.valid())
    {
        return fvc::div(phi_) + pErr_();
    }
    else
    {
        return fvc::div(phi_);
    }
}


void Foam::simpleSolver::limitField
(
    volVectorField& fld,
    const volVectorField& otherFld
)
{
    volScalarField magFld
    (
        stabilise
        (
            mag(fld),
            dimensionedScalar
            (
                "small",
                fld.dimensions(),
                SMALL
            )
        )
    );
    volScalarField clippedMagFld(min(magFld, 2*average(magFld)));
    fld *= (clippedMagFld/magFld);
}


void Foam::simpleSolver::setError
(
    simpleSolver& fineEqns,
    const meshToMesh0& fineToCoarse
)
{
    if (U0_.valid())
    {
        U0_.ref() == U_;
    }
    else
    {
        U0_ = tmp<volVectorField>(new volVectorField("U0", U_));
    }

    /*
    const volScalarField& k(mesh_.lookupObject<volScalarField>("k"));
    if (k0_.valid())
    {
        k0_() == k;
    }
    else
    {
        k0_ = tmp<volScalarField>(new volScalarField("k0", k));
    }

    const volScalarField& omega(mesh_.lookupObject<volScalarField>("omega"));
    if (omega0_.valid())
    {
        omega0_() == omega;
    }
    else
    {
        omega0_ = tmp<volScalarField>(new volScalarField("omega0", omega));
    }

    const volScalarField& nut(mesh_.lookupObject<volScalarField>("nut"));
    if (nut0_.valid())
    {
        nut0_() == nut;
    }
    else
    {
        nut0_ = tmp<volScalarField>(new volScalarField("nut0", nut));
    }
    */

    //volScalarField& nut
    //(
    //    const_cast<volScalarField&>(mesh_.lookupObject<volScalarField>("nut"))
    //);
    //nut = fineToCoarse.mapTgtToSrc
    //(fineEqns.mesh_.lookupObject<volScalarField>("nut"), plusEqOp<scalar>());

    tmp<volVectorField> thisRes(Ures(false));
    DebugVar(gMax(thisRes()));
    DebugVar(gMin(thisRes()));
    DebugVar(gAverage(thisRes()));

    tmp<volVectorField> fineRes(fineEqns.Ures());
    DebugVar(gMax(fineRes()));
    DebugVar(gMin(fineRes()));
    DebugVar(gAverage(fineRes()));

    // Note: make sure interpolation is zero where there is no overlap
    volVectorField interpolatedFineRes("interpolatedRes", thisRes());
    volScalarField weightField
    (
        IOobject
        (
            "weightField",
            fineToCoarse.toMesh().time().timeName(),
            fineToCoarse.toMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fineToCoarse.toMesh(),
        dimensionedScalar("one", dimless, 1.0)
    );
    {
        interpolatedFineRes == dimensionedVector
        (
            "zero",
            interpolatedFineRes.dimensions(),
            vector::zero
        );
        fineToCoarse.interpolate
        (
            interpolatedFineRes,
            fineRes,
            meshToMesh0::INTERPOLATE,
            plusEqOp<vector>()
        );
        forAll(fineToCoarse.inverseDistanceWeights(), coarsei)
        {
            if (fineToCoarse.inverseDistanceWeights()[coarsei].size() <= 0)
            {
                Pout<< "Coarse cell :" << coarsei
                    << " at:" << fineToCoarse.toMesh().cellCentres()[coarsei]
                    << " fine celli:"
                    << fineToCoarse.cellAddressing()[coarsei]
                    << " weights:"
                    << fineToCoarse.inverseDistanceWeights()[coarsei]
                    << endl;

                weightField[coarsei] = 0.0;
            }
        }

        // Limit field
        //limitField(interpolatedFineRes, thisRes());
    }
    DebugVar(gMax(interpolatedFineRes));
    DebugVar(gMin(interpolatedFineRes));
    DebugVar(gAverage(interpolatedFineRes));

    Uerr_ = weightField*(interpolatedFineRes - thisRes);
    Uerr_.ref().rename("Uerr");

    if (fineToCoarse.toMesh().time().outputTime())
    {
        Uerr_().write();
    }

    DebugVar(max(Uerr_()));
    DebugVar(min(Uerr_()));
    DebugVar(average(Uerr_()));


    // pErr_ =
    //     fineToCoarse.interpolate(fineEqns.pRes(), plusEqOp<scalar>())
    //   - pRes(false);
}


void Foam::simpleSolver::correct
(
    const simpleSolver& coarseEqns,
    const meshToMesh0& coarseToFine
)
{
    //volVectorField Ucorr
    //(
    //    coarseToFine.interpolate
    //    (
    //        coarseEqns.Ucorr(),
    //        meshToMesh0::INTERPOLATE,
    //        eqOp<vector>()
    //    )
    //);

    // Construct with correct bc
    volVectorField Ucorr("coarseUcorr", U_);
    {
        // Make sure that any unmapped bit gets zero correction.
        Ucorr == dimensionedVector
        (
            "zero",
            Ucorr.dimensions(),
            vector::zero
        );
//         coarseToFine.interpolate
//         (
//             Ucorr,
//             coarseEqns.Ucorr(),
//             meshToMesh0::INTERPOLATE,
//             eqOp<vector>()
//         );
        DebugVar(gMax(Ucorr()));
        DebugVar(gMin(Ucorr()));
        DebugVar(gAverage(Ucorr()));

        forAll(coarseToFine.inverseDistanceWeights(), finei)
        {
            if (coarseToFine.inverseDistanceWeights()[finei].size() <= 0)
            {
                Pout<< "Fine cell :" << finei
                    << " at:" << coarseToFine.toMesh().cellCentres()[finei]
                    << " coarse celli:"
                    << coarseToFine.cellAddressing()[finei]
                    << " weights:"
                    << coarseToFine.inverseDistanceWeights()[finei] << endl;

                Ucorr[finei] = vector::zero;
            }
        }

        // Limit field
        //limitField(Ucorr, U_);

        DebugVar(gMax(Ucorr()));
        DebugVar(gMin(Ucorr()));
        DebugVar(gAverage(Ucorr()));
        DebugVar(gMax(U_));
        DebugVar(gMin(U_));
        DebugVar(gAverage(U_));

        if (coarseToFine.toMesh().time().outputTime())
        {
            Ucorr.write();
        }
    }

    U_ += Ucorr;
    phi_ += (mesh_.Sf() & fvc::interpolate(Ucorr));

    // p_ += coarseToFine.interpolate(coarseEqns.pcorr(), plusEqOp<scalar>());

    // volScalarField& k(const_cast<volScalarField&>
    // (mesh_.lookupObject<volScalarField>("k")));
    // k += coarseToFine.interpolate(coarseEqns.kcorr(), plusEqOp<scalar>());

    // volScalarField& omega(const_cast<volScalarField&>
    ///(mesh_.lookupObject<volScalarField>("omega")));
    // omega +=
    //    coarseToFine.interpolate(coarseEqns.omegacorr(), plusEqOp<scalar>());

    // volScalarField& nut(const_cast<volScalarField&>
    // (mesh_.lookupObject<volScalarField>("nut")));
    // nut +=
    //     coarseToFine.interpolate(coarseEqns.nutcorr(), plusEqOp<scalar>());
    // nut = max(nut, dimensionedScalar("nut", nut.dimensions(), 0));
}


// ************************************************************************* //
