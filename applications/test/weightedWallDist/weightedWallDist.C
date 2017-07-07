/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    writePatchDistance

Description
    Write topological distance to wall

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "fvMesh.H"
#include "wallPointData.H"
#include "FaceCellWave.H"
#include "fvcGrad.H"
#include "localMin.H"
#include "localMax.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("patches");

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Find set of patches from the list of regular expressions provided
    const wordReList patches((IStringStream(args[1])()));
    const labelHashSet patchSet(pbm.patchSet(patches));

    if (!patchSet.size())
    {
        FatalErrorInFunction
            << "Cannot find any patches in set " << patches << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    // Load the wall distance field
    Info<< "Reading field alpha.water\n" << endl;
    volScalarField yWall
    (
        IOobject
        (
            "yWall",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // Find locations where the gradient is too low.
    volScalarField magGrad("magGrad", mag(fvc::grad(yWall)));
    magGrad.write();

    surfaceScalarField maxYWall(localMax<scalar>(mesh).interpolate(yWall));

    surfaceScalarField maxGrad
    (
        "maxGrad",
        localMax<scalar>(mesh).interpolate(magGrad)
    );
    surfaceScalarField minGrad
    (
        "minGrad",
        localMin<scalar>(mesh).interpolate(magGrad)
    );

    // Start of changes
    DynamicList<label> seedFaces(mesh.nFaces());
    DynamicList<wallPointData<scalar>> seedData(mesh.nFaces());

    // Field on cells and faces.
    List<wallPointData<scalar>> cellData(mesh.nCells());
    List<wallPointData<scalar>> faceData(mesh.nFaces());

    // Start of changes
    forAll(maxGrad, facei)
    {
        if (maxGrad[facei] > 0.5 && minGrad[facei] < 0.5)
        {
            seedFaces.append(facei);
            seedData.append
            (
                wallPointData<scalar>
                (
                    mesh.Cf()[facei],
                    maxYWall[facei],
                    0.0
                )
            );
        }
    }
    forAll(maxGrad.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pMax = maxGrad.boundaryField()[patchi];
        const fvsPatchScalarField& pMin = minGrad.boundaryField()[patchi];
        const fvsPatchScalarField& pYWall = maxYWall.boundaryField()[patchi];
        forAll(pMax, i)
        {
            if (pMax[i] > 0.5 && pMin[i] < 0.5)
            {
                label facei = pMax.patch().start()+i;
                seedFaces.append(facei);
                seedData.append
                (
                    wallPointData<scalar>
                    (
                        mesh.Cf().boundaryField()[patchi][i],
                        pYWall[i],
                        0
                    )
                );
            }
        }
    }

DebugVar(seedFaces);


    // Propagate information inwards
    FaceCellWave<wallPointData<scalar>> deltaCalc
    (
        mesh,
        seedFaces,
        seedData,
        faceData,
        cellData,
        mesh.globalData().nTotalCells()+1
    );

    volScalarField fld
    (
        IOobject
        (
            "interfaceHeight",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("large", dimless, mesh.globalData().nTotalCells()+1)
    );
    forAll(cellData, celli)
    {
        if (cellData[celli].valid(deltaCalc.data()))
        {
            fld[celli] = cellData[celli].data();
        }
    }
    forAll(fld.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& fvp = fld.boundaryFieldRef()[patchi];
        scalarField patchVals(fvp.size(), 0.0);

        forAll(patchVals, i)
        {
            label facei = fvp.patch().start()+i;
            if (faceData[facei].valid(deltaCalc.data()))
            {
                patchVals[i] = faceData[facei].data();
            }
        }

        fvp == patchVals;
    }

    Info<< "Writing " << fld.name() << " with interfaceHeight" << endl;
    fld.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
