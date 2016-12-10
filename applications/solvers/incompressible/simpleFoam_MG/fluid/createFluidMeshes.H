    const wordList fluidNames(rp["fluid"]);

    PtrList<fvMesh> fluidRegions(fluidNames.size());
    PtrList<meshToMesh0> fineToCoarseMappers(fluidNames.size());
    PtrList<meshToMesh0> coarseToFineMappers(fluidNames.size());

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
    forAll(fluidRegions, i)
    {
        if (i < fluidRegions.size()-1)
        {
            fineToCoarseMappers.set
            (
                i,
                new meshToMesh0
                (
                    fluidRegions[i],
                    fluidRegions[i+1]
                    //meshToMesh::imCellVolumeWeight
                )
            );
            coarseToFineMappers.set
            (
                i,
                new meshToMesh0
                (
                    fluidRegions[i+1],
                    fluidRegions[i]
                )
            );
        }
    }
