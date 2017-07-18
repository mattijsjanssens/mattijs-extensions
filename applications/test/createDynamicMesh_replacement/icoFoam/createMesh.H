    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    IOdictionary dict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            "",         // regionName
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );
    dict.add("dynamicFvMesh", staticFvMesh::typeName);

    autoPtr<dynamicFvMesh> meshPtr
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                dynamicFvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            ),
            dict
        )
    );

    dynamicFvMesh& mesh = meshPtr();
