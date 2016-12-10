    const wordList fluidNames(rp["fluid"]);

    PtrList<fvMesh> fluidRegions(fluidNames.size());
    PtrList<meshToMesh> fluidMappers(fluidNames.size());

    forAll(fluidNames, i)
    {
        Info<< "Create fluid mesh for region " << fluidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        fluidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    fluidNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
    for (label i = 0; i < fluidRegions.size()-1; i++)
    {
        fluidMappers.set
        (
            i,
            new meshToMesh
            (
                fluidRegions[i],
                fluidRegions[i+1],
                meshToMesh::imCellVolumeWeight
            )
        );
    }
