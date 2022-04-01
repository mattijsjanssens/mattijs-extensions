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

autoPtr<Function1<scalar>> generateShapeFunction()
{
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
    return Function1<scalar>::New("scale", dict);
}


autoPtr<mapDistribute> pointVolInterpolation
(
    const polyMesh& mesh,
    limitedDistanceData<label>::weightFunction& fun,

    const pointField& position,
    const labelList& cellID,
    const scalarField& radius,

    DynamicList<label>& collectedParticles,
    DynamicList<label>& collectedCells,
    DynamicList<scalar>& collectedWeights
)
{
    const globalIndex globalParticles(position.size());

    List<limitedDistanceData<label>> cellData(mesh.nCells());
    List<limitedDistanceData<label>> faceData(mesh.nFaces());

    // Note: initial size how?
    DynamicList<label> seedFaces(position.size());
    DynamicList<limitedDistanceData<label>> seedInfo(position.size());

    forAll(position, particlei)
    {
        const scalar overlap = fun.weight
        (
            mesh,
            cellID[particlei],
            position[particlei],
            radius[particlei]
        );

        const point& cc = mesh.cellCentres()[cellID[particlei]];

        Pout<< "For particle:" << particlei
            << " at:" << position
            << " in cell:" << cc
            << " have r:" << radius[particlei]
            << " have weight:" << overlap
            << endl;

        const limitedDistanceData<label> pInfo
        (
            scalarList(1, overlap),
            scalarList(1, radius[particlei]),
            pointList(1, position[particlei]),
            labelList(1, globalParticles.toGlobal(particlei))
        );

        const cell& cFaces = mesh.cells()[cellID[particlei]];
        for (const label facei : cFaces)
        {
            const bool changed = faceData[facei].updateFace
            (
                mesh,
                facei,
                cellID[particlei],
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

    // Collect all visited cells
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
    return autoPtr<mapDistribute>
    (
        new mapDistribute(globalParticles, collectedParticles, compactMap)
    );
}


void normaliseWeights
(
    const mapDistribute& map,
    const label nParticles,
    const List<label>& collectedParticles,
    scalarList& collectedWeights
)
{
    // Get local data in slot order
    scalarField sumWeights(nParticles, 0.0);
    forAll(collectedParticles, i)
    {
        const label particlei = collectedParticles[i];
        sumWeights[particlei] += collectedWeights[i];
    }
    // Send to originating particle and sum
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>::null(),
        nParticles,
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        sumWeights,
        scalar(0.0),
        plusEqOp<scalar>(),
        flipOp()
    );

    // Send back to cells and divide
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
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"

    // Generate a test cloud
    // ~~~~~~~~~~~~~~~~~~~~~

    passiveParticleCloud particles
    (
        mesh,
        "banana",
        IDLList<passiveParticle>()
    );
    {
        Pout<< "Adding particles." << endl;
        particles.addParticle
        (
            new passiveParticle(mesh, point(0.5001, 0.5001, 0.1001))
        );
        particles.addParticle
        (
            new passiveParticle(mesh, point(0.4001, 0.4001, 0.3001))
        );
    }


    // Extract particle data into lists
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // (radius data should really be on particle)
    scalarField radius(particles.size());
    pointField position(particles.size());
    labelList cellID(particles.size());
    {
        label particlei = 0;
        for (const auto& p : particles)
        {
            position[particlei] = p.position();
            cellID[particlei] = p.cell();
            radius[particlei] = 0.15*(particlei+1);
            particlei++;
        }
    }
    const globalIndex globalParticles(particles.size());


    // Generate weight function
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Currently : approx overlap volume
    auto scalePtr(generateShapeFunction());
    limitedDistanceData<label>::weightFunction fun(scalePtr());


    // Get list of cells+weights+originating particle slot
    DynamicList<label> spreadParticles;
    DynamicList<label> spreadCells;
    DynamicList<scalar> spreadWeights;
    autoPtr<mapDistribute> mapPtr
    (
        pointVolInterpolation
        (
            mesh,
            fun,

            position,
            cellID,
            radius,

            spreadParticles,
            spreadCells,
            spreadWeights
        )
    );
    const auto& map = mapPtr();


    if (true)
    {
        // Make sure all weights add up to 1. Build into pointVolInterpolation
        // routine?
        normaliseWeights
        (
            map,
            particles.size(),
            spreadParticles,
            spreadWeights
        );
    }


    // Test: pull data from particles to spread
    {
        map.distribute(position);
        forAll(spreadParticles, i)
        {
            const label celli = spreadCells[i];
            const label sloti = spreadParticles[i];
            Pout<< "At cell:" << celli << " at:" << mesh.cellCentres()[celli]
                << " have origin:" << position[sloti]
                << " with weight:" << spreadWeights[i]
                << endl;
        }
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
    forAll(spreadParticles, i)
    {
        const label celli = spreadCells[i];
        vsf[celli] += spreadWeights[i];
    }
    vsf.correctBoundaryConditions();

    Info<< "Writing weight volScalarField to " << runTime.value()
        << endl;

    vsf.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
