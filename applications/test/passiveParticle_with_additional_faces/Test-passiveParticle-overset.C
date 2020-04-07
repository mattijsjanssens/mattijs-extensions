/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 M. Janssens
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
    Test 'interpolating' particles across overset interpolation. This is
    a cell to cell interpolation. Does:
    - collect particles per cell (see e.g. DSMCCloud::cellOccupancy()).
    - send ('interpolate') thier locations to overlapping cells (from donor
      cells to acceptor cells in overset speak).
    - find locally. Send back -1 or global cell number if found.
    - actually send over selected particles (remove from local clould) and merge
      received particles into local cloud
    The alternative is to do it in one step: send over the particles the to do
    the testing. However each acceptor might have say 7 donor cells
    (central donor + 6 neighbours)
    so there potentially a lot of data to transfer just to find who contains
    the particle.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "passiveParticleCloud.H"
#include "OBJstream.H"
#include "cellCellStencil.H"
#include "cellCellStencilObject.H"
#include "dynamicOversetFvMesh.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void collectParticlePositions
(
    const polyMesh& mesh,
    const passiveParticleCloud& particles,
    const bitSet& isDonor,  // cells to potentially send over particles
    List<pointList>& cellPositions
)
{
    // Build cell occupancy (particle position only)

    // Collect particle locations
    labelList nParticles(mesh.nCells(), 0);
    for (const passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (isDonor[celli])
        {
            nParticles[celli]++;
        }
    }

    cellPositions.setSize(mesh.nCells());

    for (const passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (isDonor[celli])
        {
            cellPositions[celli].setSize(nParticles[celli]);
        }
    }
    nParticles = 0;
    for (const passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (isDonor[celli])
        {
            label& slot = nParticles[celli];
            cellPositions[celli][slot] = p.position();
            slot++;
        }
    }
}


void findPositions
(
    const polyMesh& mesh,
    const globalIndex& globalCells,
    const mapDistribute& map,
    List<pointList>& cellPositions,
    const labelList& cellIDs,       // acceptor cells
    const labelListList& stencil,   // donor cells/slots per acceptor

    labelListList& cellDestinations // per donor cell, per particle the acceptor
)
{
    // Collect remote contributions. Now each acceptor can access
    // through the stencil the particle positions in each of its
    // donors.
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        mesh.nCells(),
        map.subMap(),
        false,
        map.constructMap(),
        false,
        cellPositions,
        pointList(),                            // nullValue
        ListOps::appendEqOp<point>(),
        flipOp(),                               // negateOp
        UPstream::msgType(),
        map.comm()
    );


    cellDestinations.setSize(map.constructSize());

    // See what we can find locally. Mark unfound slots
    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        const labelList& nbrs = stencil[celli];

        for (const label sloti : nbrs)
        {
            const pointList& particles = cellPositions[sloti];
            labelList& destinations = cellDestinations[sloti];
            if (particles.size())
            {
                // Size cell array (if not yet sized). Initialise to -1 to
                // indicate not yet found.
                destinations.setSize(particles.size(), -1);
                forAll(particles, parti)
                {
                    const point& pt = particles[parti];
                    // Check that slot(=particle) hasn't been found yet
                    if (destinations[parti] == -1)
                    {
                        if (mesh.pointInCell(pt, celli))
                        {
                            Pout<< "    Found " << pt
                                << " in cell:" << celli
                                << " at:" << mesh.cellCentres()[celli]
                                << endl;
                            destinations[parti] = globalCells.toGlobal(celli);
                        }
                    }
                }
            }
        }
    }

    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        mesh.nCells(),
        map.constructMap(),
        false,
        map.subMap(),
        false,
        cellDestinations,
        labelList(),                            // nullValue
        ListOps::appendEqOp<label>(),
        flipOp(),                               // negateOp
        UPstream::msgType(),
        map.comm()
    );
}


void extractTransferParticles
(
    const polyMesh& mesh,
    const globalIndex& globalCells,
    const bitSet& isDonor,
    const labelListList& cellDestinations,
    passiveParticleCloud& particles,
    List<IDLList<passiveParticle>>& particleTransferLists,
    List<DynamicList<point>>& positionTransferLists
)
{    
    // Allocate transfer buffers
    particleTransferLists.setSize(Pstream::nProcs());
    positionTransferLists.setSize(Pstream::nProcs());

    // Loop over particles in same order
    labelList nParticles(mesh.nCells(), 0);
    for (passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (isDonor[celli] && cellDestinations[celli].size())
        {
            label& slot = nParticles[celli];
            const label acceptorCelli = cellDestinations[celli][slot];
            slot++;

            if (acceptorCelli != -1)
            {
                const label proci = globalCells.whichProcID(acceptorCelli);
                const label remoteCelli =
                    globalCells.toLocal(proci, acceptorCelli);
                const point pos(p.position());
                //Pout<< "** Transferring particle " << pos
                //    << " from cell:" << celli
                //    << " to proc:" << proci
                //    << " remote cell:" << remoteCelli << endl;
                passiveParticle* pPtr = particles.remove(&p);
                pPtr->cell() = remoteCelli;
                particleTransferLists[proci].append(pPtr);
                positionTransferLists[proci].append(pos);
            }
        }
    }
}


void distributeParticles
(
    const polyMesh& mesh,
    const List<IDLList<passiveParticle>>& particleTransferLists,
    const List<DynamicList<point>>& positionTransferLists,
    passiveParticleCloud& particles
)
{
    // Allocate transfer buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    forAll(particleTransferLists, proci)
    {
        if (particleTransferLists[proci].size())
        {
            UOPstream particleStream(proci, pBufs);
            particleStream
                << positionTransferLists[proci]
                << particleTransferLists[proci];
        }
    }
    labelList allNTrans(Pstream::nProcs());
    pBufs.finishedSends(allNTrans);

    forAll(allNTrans, neighbProci)
    {
        label nRec = allNTrans[neighbProci];

        if (nRec)
        {
            UIPstream particleStream(neighbProci, pBufs);

            Pout<< "From processor " << neighbProci
                << " receiving " << nRec << " bytes" << endl;

            const pointList newPos(particleStream);
            IDLList<passiveParticle> newParticles
            (
                particleStream,
                typename passiveParticle::iNew(mesh)
            );

            label pI = 0;

            for (passiveParticle& newp : newParticles)
            {
                // Particle already has correct cell but incorrect position
                newp.relocate(newPos[pI], newp.cell());
                particles.append(newParticles.remove(&newp));
                pI++;
            }
        }
    }
}


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

    const pointField& cellCentres = mesh.cellCentres();

    // Seeding random particles
    Random rndGen(0);
    for (label i = 0; i < 1000; i++)
    {
        const label celli = rndGen.position(0, mesh.nCells()-1);
        const cell& cFaces = mesh.cells()[celli];
        const boundBox bb(cFaces.points(mesh.faces(), mesh.points()), false);

        const point pt(rndGen.position(bb.min(), bb.max()));
        particles.addParticle
        (
            new passiveParticle(mesh, pt, celli)
        );
    }

    // Check
    for (const passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (!mesh.pointInCell(p.position(), celli))
        {
            const cell& c = mesh.cells()[celli];
            const boundBox bb(c.points(mesh.faces(), mesh.points()), false);
            FatalErrorInFunction << "Particle at:" << p.position()
                << " is not in cell " << celli
                << " cc:" << cellCentres[celli]
                << " bb:" << bb << exit(FatalError);
        }
    }
    Pout<< "Particles:" << particles.size() << endl;

    runTime.printExecutionTime(Info);

    //++runTime;
    //Pout<< "Writing particles to time " << runTime.timeName() << endl;
    //particles.write();


    const globalIndex globalCells(mesh.nCells());

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


    // Collect particle locations (for selected cells)
    List<pointList> cellPositions;
    collectParticlePositions
    (
        mesh,
        particles,
        isDonor,        // cells to collect particles for
        cellPositions
    );


    // Find per acceptor cell (in cellIDs), per position, the global
    // cell containing it. Note: extends cellPositions with remote positions.
    labelListList cellDestinations;
    findPositions
    (
        mesh,
        globalCells,
        map,
        cellPositions,
        cellIDs,            // acceptor cells
        stencil,            // donor cells/slots per acceptor

        cellDestinations    // per donor cell, per particle the acceptor
    );


    for (label celli = 0; celli < mesh.nCells(); celli++)
    {
        const pointList& cellSamples = cellPositions[celli];
        const labelList& acceptors = cellDestinations[celli];
        if (acceptors.size())
        {
            Pout<< "For cell:" << celli << " at:" << cellCentres[celli]
                << " have samples:" << cellSamples
                << " acceptors:" << acceptors << endl;
        }
    }


    // Split particles into
    // - remaining (in original Cloud)
    // - to be sent to other processor
    // For the sent particles we need to store both the particle itself
    // and its original location since the particle only stores the
    // barycentric coordinates (which are no longer relevant in the new
    // cell)

    // Allocate transfer buffers
    List<IDLList<passiveParticle>> particleTransferLists(Pstream::nProcs());
    List<DynamicList<point>> positionTransferLists(Pstream::nProcs());
    // Filter particles into local ones and ones that need to transfer    
    extractTransferParticles
    (
        mesh,
        globalCells,
        isDonor,
        cellDestinations,

        particles,
        particleTransferLists,
        positionTransferLists
    );

    Pout<< "** REMAINNG Particles:" << particles.size() << endl;

    // Send over transferlists and append to particles. Make sure to use
    // old location to update barycentric coordinates to the new cell
    distributeParticles
    (
        mesh,
        particleTransferLists,
        positionTransferLists,
        particles
    );


    Pout<< "** All Particles:" << particles.size() << endl;

    // Check
    for (const passiveParticle& p : particles)
    {
        const label celli = p.cell();
        if (!mesh.pointInCell(p.position(), celli))
        {
            const cell& c = mesh.cells()[celli];
            const boundBox bb(c.points(mesh.faces(), mesh.points()), false);
            FatalErrorInFunction << "Particle at:" << p.position()
                << " is not in cell " << celli
                << " cc:" << cellCentres[celli]
                << " bb:" << bb << exit(FatalError);
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
