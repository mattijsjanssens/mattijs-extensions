    // Momentum predictor

    #include "assembleUEqn.H"

    if (simple.momentumPredictor())
    {
        tmp<fvVectorMatrix> tUpEqn(UEqn == -fvc::grad(p));

        fvVectorMatrix& UpEqn = tUpEqn.ref();

        UpEqn += UDefect;

        solve(UpEqn);

        fvOptions.correct(U);
    }
