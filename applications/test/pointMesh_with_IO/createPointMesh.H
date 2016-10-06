Foam::Info
    << "Create pointMesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

Foam::pointMesh pMesh
(
    Foam::IOobject
    (
        mesh.name(),                // polyMesh region
        mesh.facesInstance(),       // polyMesh topology instance
        runTime,
        Foam::IOobject::MUST_READ
    ),
    mesh
);
