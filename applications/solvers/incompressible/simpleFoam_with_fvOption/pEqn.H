{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF);

    {
        volScalarField divPhiHbyA(fvc::div(phiHbyA));
        // Non-orthogonal pressure corrector loop
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAtU(), p) == divPhiHbyA
            );

            pEqn.setReference(pRefCell, pRefValue);

            fvOptions.constrain(pEqn);

            pEqn.solve();

            if (simple.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();

                // Adjust flux on cell
                {
                    const label celli = 1001;

                    scalar massIn(0.0);
                    scalar massOut(0.0);
                    label nFaces = 0;
                    const cell& cFaces = mesh.cells()[1001];
                    for (const label facei : cFaces)
                    {
                        if (mesh.isInternalFace(facei))
                        {
                            const label own = mesh.faceOwner()[facei];
                            const bool outFlow =
                            (
                                (own == celli)
                             == (phi[facei] > 0)
                            );

                            Pout<< "    across face:" << facei
                                << " at:" << mesh.faceCentres()[facei]
                                << " have flux:" << phi[facei]
                                << " in direction:"
                                << mesh.faces()[facei].unitNormal(mesh.points())
                                << endl;

                            if (outFlow)
                            {
                                massOut += mag(phi[facei]);
                            }
                            else
                            {
                                massIn += mag(phi[facei]);
                            }
                            nFaces++;
                        }
                    }

                    DebugVar(massOut);
                    DebugVar(massIn);


                    // Distribute any difference - 50% to inflow,
                    // 50% to outflow
                    const scalar massLoss = massOut-massIn;
                    const scalar s = 0.5*massLoss/nFaces;

                    massIn = 0.0;
                    massOut = 0.0;
                    for (const label facei : cFaces)
                    {
                        if (mesh.isInternalFace(facei))
                        {
                            const label own = mesh.faceOwner()[facei];
                            const bool outFlow =
                            (
                                (own == celli)
                             == (phi[facei] > 0)
                            );

                            Pout<< "    across face:" << facei
                                << " at:" << mesh.faceCentres()[facei]
                                << " have flux:" << phi[facei]
                                << " in direction:"
                                << mesh.faces()[facei].unitNormal(mesh.points())
                                << endl;

                            if (outFlow)
                            {
                                phi[facei] -= s;
                                massOut += mag(phi[facei]);
                            }
                            else
                            {
                                phi[facei] += s;
                                massIn += mag(phi[facei]);
                            }
                        }
                    }
                    DebugVar(massOut);
                    DebugVar(massIn);

                }
            }
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAtU()*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
