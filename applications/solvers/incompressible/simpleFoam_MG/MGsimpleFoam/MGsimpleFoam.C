/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    MGsimpleFoam

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "meshToMesh.H"
#include "simpleSolver.H"
#include "PtrList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    PtrList<Time> timeLevels(1);
    PtrList<fvMesh> meshLevels(1);

    timeLevels.set
    (
        0,
        new Time(Time::controlDictName, args)
    );

    meshLevels.set
    (
        0,
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                timeLevels[0].timeName(),
                timeLevels[0],
                IOobject::MUST_READ
            )
        )
    );

    const dictionary& MGdict(meshLevels[0].fvSolution::subDict("MG"));
    const label nLevels(readLabel(MGdict.lookup("nLevels")));

DebugVar(nLevels);

    const label nFinestIter(readLabel(MGdict.lookup("nFinestIter")));
    const label nCoarsestIter(readLabel(MGdict.lookup("nCoarsestIter")));
    const label nIter(readLabel(MGdict.lookup("nIter")));

    timeLevels.resize(nLevels);
    meshLevels.resize(nLevels);

    PtrList<meshToMesh> mapLevels(nLevels - 1);
    PtrList<simpleSolver> eqnLevels(nLevels);

    for (label level=1; level<nLevels; level++)
    {
        timeLevels.set
        (
            level,
            new Time
            (
                Time::controlDictName,
                args.path()/"MG",
                "level" + name(level)
            )
        );

        meshLevels.set
        (
            level,
            new fvMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    timeLevels[level].timeName(),
                    timeLevels[level],
                    IOobject::MUST_READ
                )
            )
        );
    }

    for (label level=0; level<nLevels-1; level++)
    {
        mapLevels.set
        (
            level,
            new meshToMesh
            (
                meshLevels[level + 1],
                meshLevels[level],
                meshToMesh::imCellVolumeWeight
                //meshToMesh::imDirect
                //meshToMesh::imMapNearest
            )
        );
    }

    for (label level=0; level<nLevels; level++)
    {
        eqnLevels.set
        (
            level,
            new simpleSolver(meshLevels[level])
        );
    }

    label cycle = 0;
    while (eqnLevels[0].loop())
    {
        Info<< "V-cycle " << ++cycle << endl;

        // Restrict residuals from fine to coarse meshes
        for (label level=0; level<nLevels-1; level++)
        {
            eqnLevels[level+1].setError
            (
                eqnLevels[level],
                mapLevels[level]
            );
        }


        // Solve the coarsest mesh
        timeLevels[nLevels-1]++;
        eqnLevels[nLevels-1].solve(nCoarsestIter);

        // Prolong correction from coarse meshes to finer and solve
        for (label level=nLevels-2; level>0; level--)
        {
            Info<< "Correcting " << level << endl;
            eqnLevels[level].correct
            (
                eqnLevels[level+1],
                mapLevels[level]
            );

            timeLevels[level]++;
            eqnLevels[level].solve(nIter);
        }

        // Prolong correction to finest mesh
        eqnLevels[0].correct
        (
            eqnLevels[1],
            mapLevels[0]
        );

        // Solve the finest mesh
        timeLevels[0]++;
        eqnLevels[0].solve(nFinestIter);

        forAll(timeLevels, level)
        {
            timeLevels[level].write();
        }

        Info<< "ExecutionTime = " << timeLevels[0].elapsedCpuTime() << " s"
            << "  ClockTime = " << timeLevels[0].elapsedClockTime() << " s"
            << nl << endl;

        Info<< endl;
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
