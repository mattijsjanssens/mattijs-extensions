/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    testPassiveParticle

Description
    Test cloud of passive particles.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "passiveParticleCloud.H"
#include "OBJstream.H"
#include "cellCellStencil.H"
#include "cellCellStencilObject.H"
#include "dynamicOversetFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noFunctionObjects();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    const word cloudName = "myFirstCloud";

    // Start with empty cloud
    passiveParticleCloud particles
    (
        mesh,
        cloudName,
        IDLList<passiveParticle>()
    );
    Pout<< "Starting particles:" << particles.size() << endl;

    Pout<< "Adding particles." << endl;
    const pointField& cellCentres = mesh.cellCentres();
    forAll(cellCentres, celli)
    {
        particles.addParticle
        (
            new passiveParticle(mesh, cellCentres[celli], celli)
        );
    }

    //for (const passiveParticle& p : particles)
    //{
    //    Pout<< "    " << p.position() << " cell:" << p.cell()
    //        << " origProc:" << p.origProc()
    //        << " origId:" << p.origId()
    //        << " coord:" << p.coordinates()
    //        << " transform:" << p.currentTetTransform()
    //        << endl;
    //}

    runTime.printExecutionTime(Info);

    //++runTime;
    //Pout<< "Writing particles to time " << runTime.timeName() << endl;
    //particles.write();


    const cellCellStencil& overlap = Stencil::New(mesh);
    const labelListList& stencil = overlap.cellStencil();
    const mapDistribute& map = overlap.cellInterpolationMap();
    const labelList& cellIDs = overlap.interpolationCells();


    // Mark donor cells
    bitSet isDonor(mesh.nCells());
    for (const labelList& sendCells : map.subMap())
    {
        isDonor.set(sendCells);
    }
    //DebugVar(isDonor);

    // Collect particle locations
    labelList nParticles(mesh.nCells(), 0);
    for (const passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (isDonor[celli])
        {
            nParticles[celli]++;
            //Pout<< "Donor:" << celli << " at:" << cellCentres[celli]
            //    << " now " << nParticles[celli] << endl;
        }
    }
    List<pointList> cellOccupancy(mesh.nCells());
    for (const passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (isDonor[celli])
        {
            cellOccupancy[celli].setSize(nParticles[celli]);
        }
    }
    nParticles = 0;
    for (const passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (isDonor[celli])
        {
            cellOccupancy[celli][nParticles[celli]++] = p.position();
        }
    }
    DebugVar(cellOccupancy);

    // Collect remote contributions
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        mesh.nCells(),
        map.subMap(),
        false,
        map.constructMap(),
        false,
        cellOccupancy,
        pointList(),                            // nullValue
        ListOps::appendEqOp<point>(),
        flipOp(),                               // negateOp
        UPstream::msgType(),
        map.comm()
    );

    // See what we can find locally. Mark unfound slots
    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        const labelList& nbrs = stencil[celli];

        if (nbrs.size())
        {
            Pout<< "Acceptor cell:" << celli << " at:"
                << mesh.cellCentres()[celli] << nl;
            for (const label sloti : nbrs)
            {
                pointList& particles = cellOccupancy[sloti];
                if (particles.size())
                {
                    Pout<< "    from donor:" << sloti
                        << " potential donor particles:"
                        << particles << endl;

                    for (point& pt : particles)
                    {
                        if (pt != point::max)
                        {
                            if (!mesh.pointInCell(pt, celli))
                            {
                                Pout<< "    did not find " << pt
                                    << " in cell:" << mesh.cellCentres()[celli]
                                    << endl;
                                pt = point::max;
                            }
                            else
                            {
                                Pout<< "    Find " << pt << " in cell:"
                                    << mesh.cellCentres()[celli] << endl;
                            }
                        }
                    }
                }
            }
        }
    }

    // Send back
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        mesh.nCells(),
        map.constructMap(),
        false,
        map.subMap(),
        false,
        cellOccupancy,
        pointList(),                            // nullValue
        ListOps::appendEqOp<point>(),
        flipOp(),                               // negateOp
        UPstream::msgType(),
        map.comm()
    );

    for (label celli = 0; celli < mesh.nCells(); celli++)
    {
        const pointList& cellSamples = cellOccupancy[celli];
        if (cellSamples.size())
        {
            Pout<< "For cell:" << celli << " at:" << cellCentres[celli]
                << " have samples:" << cellSamples << endl;
        }
    }
    


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
