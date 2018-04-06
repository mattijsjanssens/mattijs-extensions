Foam::Info
    << "Create unallocated mesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::autoPtr<Foam::unallocatedFvMesh> meshPtr
(
    Foam::unallocatedFvMeshTools::newMesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    )
);
Foam::unallocatedFvMesh& mesh = meshPtr();
