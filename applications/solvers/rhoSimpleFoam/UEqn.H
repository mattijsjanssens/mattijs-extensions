    // Solve the Momentum equation

    Pout<< "BEFORE U: min:" << gMin(U) << " max:" << gMax(U) << endl;
    Pout<< "BEFORE rho: min:" << gMin(rho) << " max:" << gMax(rho) << endl;

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
    Pout<< "AFTER U: min:" << gMin(U) << " max:" << gMax(U) << endl;

