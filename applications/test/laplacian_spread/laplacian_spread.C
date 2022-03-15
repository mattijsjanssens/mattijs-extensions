/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022 M. Janssens
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
    laplacian_spread

Description
    Lagrangian-to-VOF (i.e. spreading a point function)
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "fvMesh.H"
#include "FaceCellWave.H"
#include "limitedDistanceData.H"
#include "EdgeMap.H"
#include "passiveParticleCloud.H"
#include "columnFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"

    //const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    // Find set of patches from the list of regular expressions provided
    //const word patchName((IStringStream(args[1])()));
    //label patchi = pbm.findPatchID(patchName);

    // Generate our own cloud
    passiveParticleCloud particles
    (
        mesh,
        "banana",
        IDLList<passiveParticle>()
    );
    Pout<< "Adding particles." << endl;
    particles.addParticle
    (
        new passiveParticle(mesh, point(0.5001, 0.5001, 0.1001))
    );
    particles.addParticle
    (
        new passiveParticle(mesh, point(0.4001, 0.4001, 0.3001))
    );

    // Generate indices
    List<const passiveParticle*> particleList(particles.size());
    scalarField radius(particles.size());
    
    label particlei = 0;
    for (const auto& p : particles)
    {
        particleList[particlei] = &p;
        radius[particlei] = 0.15*(particlei+1);
        particlei++;
    }

    const globalIndex globalParticles(particleList.size());

    // Generate 
    List<limitedDistanceData<label>> cellData(mesh.nCells());
    List<limitedDistanceData<label>> faceData(mesh.nFaces());

    // Shape of sphere. Needs to be zero outside of sphere.
    dictionary scaleDict;
    List<Tuple2<scalar, scalar>> vals(3);
    vals[0] = Tuple2<scalar, scalar>(0.0, 1.0);
    vals[1] = Tuple2<scalar, scalar>(1, 1.0);
    // Weird: uses streaming so truncates to writePrecision ...
    vals[2] = Tuple2<scalar, scalar>(1.00001, 0.0);
    scaleDict.add("values", vals);
    scaleDict.add("type", "table");
    dictionary dict;
    dict.add("scale", scaleDict);
    autoPtr<Function1<scalar>> scalePtr(Function1<scalar>::New("scale", dict));

    DebugVar(scalePtr().value(0.1));
    DebugVar(scalePtr().value(1));
    DebugVar(scalePtr().value(1.1));
    DebugVar(scalePtr().value(200));


    limitedDistanceData<label>::weightFunction fun(scalePtr());


    DynamicList<label> seedFaces;
    DynamicList<limitedDistanceData<label>> seedInfo;

    forAll(particleList, particlei)
    {
        const auto& p = *particleList[particlei];
        const label celli = p.cell();
        const vector position = p.position();

        const scalar overlap = fun.weight
        (
            mesh,
            celli,
            position,
            radius[particlei]
        );

        const point& cc = mesh.cellCentres()[celli];

        Pout<< "For particle:" << particlei
            << " at:" << position
            << " in cell:" << cc
            << " have r:" << radius[particlei]
            << " have weight:" << overlap
            << " have scale:"
            << scalePtr().value(mag(cc-position))
            << endl;

        const limitedDistanceData<label> pInfo
        (
            scalarList(1, overlap),
            scalarList(1, radius[particlei]),
            pointList(1, position),
            labelList(1, globalParticles.toGlobal(particlei))
        );

        const cell& cFaces = mesh.cells()[celli];
        for (const label facei : cFaces)
        {
            const bool changed = faceData[facei].updateFace
            (
                mesh,
                facei,
                celli,
                pInfo,
                SMALL,
                fun
            );
            if (changed)
            {
                seedFaces.append(facei);
                seedInfo.append(pInfo);
            }
        }
    }


    // Propagate information inwards
    FaceCellWave
    <
        limitedDistanceData<label>,
        limitedDistanceData<label>::weightFunction
    > deltaCalc
    (
        mesh,
        seedFaces,
        seedInfo,
        faceData,
        cellData,
        mesh.globalData().nTotalCells()+1,
        fun
    );


    // - make untemplated - data is always global particle index
    // - better volume overlap calculation
    // - allow normalisation bypass for
    //      -step scaling
    //      -exact overlap calculation


    // Normalisation of weights


    // Collect all visited cells
    DynamicList<label> collectedParticles;
    DynamicList<scalar> collectedWeights;
    DynamicList<label> collectedCells;

    forAll(cellData, celli)
    {
        const auto& data = cellData[celli].data();
        const auto& weights = cellData[celli].weights();
        forAll(data, i)
        {
            collectedParticles.append(data[i]); // global index of particle
            collectedWeights.append(weights[i]);
            collectedCells.append(celli);
        }
    }

    Pout<< "Have collected cells:" << collectedCells.size()
        << endl;

    collectedParticles.shrink();
    collectedWeights.shrink();
    collectedCells.shrink();

    List<Map<label>> compactMap;
    mapDistribute map(globalParticles, collectedParticles, compactMap);

    //// Test: pull data from particles to collected
    //{
    //    pointField position(particleList.size());
    //    forAll(particleList, particlei)
    //    {
    //        position[particlei] = particleList[particlei]->position();
    //    }
    //    map.distribute(position);
    //    forAll(collectedParticles, i)
    //    {
    //        const label celli = collectedCells[i];
    //        const label particlei = collectedParticles[i];
    //        Pout<< "At cell:" << celli << " at:" << mesh.cellCentres()[celli]
    //            << " have origin:" << position[particlei]
    //            << " with weight:" << collectedWeights[i]
    //            << endl;
    //    }
    //}

    // Sum weights onto particle
    scalarField sumWeights(particles.size(), 0.0);
    forAll(collectedParticles, i)
    {
        const label particlei = collectedParticles[i];
        sumWeights[particlei] += collectedWeights[i];
    }
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>::null(),
        particles.size(),           // constructSize
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        sumWeights,
        scalar(0.0),
        plusEqOp<scalar>(),
        flipOp()
    );

    // Send back and divide
    map.distribute(sumWeights);
    forAll(collectedParticles, i)
    {
        //const label celli = collectedCells[i];
        const label particlei = collectedParticles[i];
        //Pout<< "At cell:" << celli << " at:" << mesh.cellCentres()[celli]
        //    << " have origin:" << position[particlei]
        //    << " with weight:" << collectedWeights[i]
        //    << " sum:" << sumWeights[particlei]
        //    << endl;
       collectedWeights[i] /= sumWeights[particlei];
    }

    // Extract as volScalarField
    volScalarField vsf
    (
        IOobject
        (
            "weight",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("weight", dimLength, 0.0)
    );
    forAll(collectedParticles, i)
    {
        const label celli = collectedCells[i];
        vsf[celli] += collectedWeights[i];
    }
    vsf.correctBoundaryConditions();

    Info<< "Writing weight volScalarField to " << runTime.value()
        << endl;

    vsf.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
