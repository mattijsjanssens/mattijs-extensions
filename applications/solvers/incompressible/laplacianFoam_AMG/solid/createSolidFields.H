    // Initialise solid field pointer lists
    PtrList<volScalarField> Ts(solidRegions.size());
    PtrList<fv::options> solidHeatSources(solidRegions.size());
    PtrList<dimensionedScalar> DTs(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to T\n" << endl;
        Ts.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding fvOptions\n" << endl;
        solidHeatSources.set
        (
            i,
            new fv::options(solidRegions[i])
        );

        Info<< "Reading transportProperties\n" << endl;

        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );


        Info<< "Reading diffusivity DT\n" << endl;

        DTs.set
        (
            i,
            new dimensionedScalar
            (
                transportProperties.lookup("DT")
            )
        );
    }
