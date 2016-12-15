// Initialise fluid field pointer lists
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> pFluidRes(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<volVectorField> UFluidRes(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());

PtrList<singlePhaseTransportModel> laminarTransport(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulence(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size(), 0);
List<scalar> pRefValueFluid(fluidRegions.size(), 0.0);

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to p\n" << endl;
    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );
    pFluidRes.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("zero", pFluid[i].dimensions()/dimTime, 0.0)
        )
    );

    Info<< "    Adding to U\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );
    UFluidRes.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedVector("zero", UFluid[i].dimensions()/dimTime, vector::zero)
        )
    );

    Info<< "    Adding to phi\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(UFluid[i])
        )
    );

    Info<< "    Adding to transport\n" << endl;
    laminarTransport.set
    (
        i,
        new singlePhaseTransportModel
        (
            UFluid[i],
            phiFluid[i]
        )
    );


    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            UFluid[i],
            phiFluid[i],
            laminarTransport[i]
        ).ptr()
    );

    const dictionary& simpleDict =
        fluidRegions[i].solutionDict().subDict("SIMPLE");

    setRefCell
    (
        pFluid[i],
        simpleDict,
        pRefCellFluid[i],
        pRefValueFluid[i]
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    turbulence[i].validate();
}
