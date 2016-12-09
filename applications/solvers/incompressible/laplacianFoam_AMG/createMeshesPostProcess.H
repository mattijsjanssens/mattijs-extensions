#include "createMeshes.H"

if (!solidRegions.size())
{
    FatalErrorIn(args.executable())
        << "No solid meshes present" << exit(FatalError);
}

fvMesh& mesh = solidRegions[0];
