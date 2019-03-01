rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
tUEqn.clear();

bool closedVolume = false;

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

volScalarField rhorAtU("rhorAtU", rho*rAtU);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAtU, MRF);


Pout<< "before p: min:" << min(p) << " max:" << max(p) << endl;
Pout<< "before rho: min:" << min(rho) << " max:" << max(rho) << endl;
Pout<< "before U: min:" << min(U) << " max:" << max(U) << endl;



if (simple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );

    phiHbyA +=
        fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf()
      - fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

    HbyA -= (rAU - rAtU)*fvc::grad(p);

    const volScalarField oldP("oldP", p);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAtU, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        // Relax the pressure equation to maintain diagonal dominance
        pEqn.relax();

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            // underrelax
            Info<< "p min:" << min(p) << " oldP min:" << min(oldP) << endl;
            //const scalar pMinFactor(min(p).value()/min(oldP).value());
            //Info<< "Detected p change:" << pMinFactor << endl;
            //if (pMinFactor < 0.1)
            //{
            //    if (pMinFactor > 0.0)
            //    {
            //        scalar alpha = pMinFactor/0.1;
            //        Info<< "*** Underrleax:" << alpha << endl;
            //        p = alpha*p+(1.0-alpha)*oldP;
            //        phi = phiHbyA + alpha*pEqn.flux();
            //    }
            //    else
            //    {
            //        scalar alpha = 0.1;
            //        Info<< "*** Underrleax:" << alpha << endl;
            //        p = alpha*p+(1.0-alpha)*oldP;
            //        phi = phiHbyA + alpha*pEqn.flux();
            //    }
            //}
            //else
            //{
            //    phi = phiHbyA + pEqn.flux();
            //}
            label minCelli = findMin(p);
            Pout<< "minCelli:" << minCelli << " p:" << p[minCelli]
                << " oldP:" << oldP[minCelli] << endl;
            if (p[minCelli] < 1e4)
            {
                //const scalar alpha =
                //    (oldP[minCelli]-1e4)
                //   /(oldP[minCelli]-p[minCelli]);
                const scalar alpha =
                    (min(oldP).value()-1e4)
                   /(min(oldP).value()-min(p).value());

                p = alpha*p+(1.0-alpha)*oldP;
                phi = phiHbyA + alpha*pEqn.flux();

                Info<< "*** Underrleax:" << alpha
                    << " afterards p:" << p[minCelli]
                    << " addedflux:" << max(mag(alpha*pEqn.flux())) << endl;
            }
            else
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
}
else
{
    closedVolume = adjustPhi(phiHbyA, U, p);

    phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU)*fvc::grad(p);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rhorAtU, p)
          ==
            fvOptions(psi, p, rho.name())
        );

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

// The incompressibe form of the continuity error check is appropriate for
// steady-state compressible also.
#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

bool pLimited = pressureControl.limit(p);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

if (pLimited || closedVolume)
{
    p.correctBoundaryConditions();
}

rho = thermo.rho();

if (!simple.transonic())
{
    rho.relax();
}


Pout<< "p: min:" << min(p) << " max:" << max(p) << endl;
Pout<< "rho: min:" << min(rho) << " max:" << max(rho) << endl;
Pout<< "U: min:" << min(U) << " max:" << max(U) << endl;

