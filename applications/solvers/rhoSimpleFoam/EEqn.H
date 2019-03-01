{
    volScalarField& he = thermo.he();

    const volScalarField oldHe("oldHe", he);

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

{
    Pout<< "BEFORE he: min:" << min(oldHe) << " max:" << max(oldHe) << endl;
    Pout<< "he: min:" << min(he) << " max:" << max(he) << endl;
    const scalar heMinFactor(min(he).value()/min(oldHe).value());
    Info<< "Detected he change:" << heMinFactor << endl;

    label minCelli = findMin(he);
    Pout<< "minCelli:" << minCelli << " he:" << he[minCelli]
        << " oldHe:" << oldHe[minCelli] << endl;
    if (he[minCelli] < 1e5)
    {
        //const scalar alpha =
        //    (oldHe[minCelli]-1e5)
        //   /(oldHe[minCelli]-he[minCelli]);

        const scalar alpha =
            (min(oldHe).value()-1e5)
           /(min(oldHe).value()-min(he).value());


        Info<< "**** he alpha:" << alpha << endl;
        he = alpha*he + (1.0-alpha)*oldHe;
        Pout<< "CORRECTED he: min:" << min(he) << " max:" << max(he) << endl;
    }
}

    fvOptions.correct(he);

    thermo.correct();
}
