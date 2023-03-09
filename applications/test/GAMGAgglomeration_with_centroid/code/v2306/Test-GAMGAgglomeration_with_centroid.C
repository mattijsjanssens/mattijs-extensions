/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 M. Janssens
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
    Test-GAMGAgglomeration

Description
    Test application for GAMG agglomeration. Hardcoded to expect GAMG on p.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "GAMGAgglomeration.H"
#include "OBJstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool closer(const point& at, const point& origin, point& nearest)
{
    if (magSqr(at-origin) < magSqr(at-nearest))
    {
        // Origin is closer
        nearest = origin;
        return true;
    }
    return false;
}


void meshWaveNearest
(
    const lduMesh& mesh,
    const pointField& faceCentres,
    const pointField& cellCentres,
    const labelUList& seedFaces,
    const List<point>& seedInfo,
    List<point>& allFaceInfo,
    List<point>& allCellInfo
)
{
    const lduAddressing& addr = mesh.lduAddr();
    const label* const __restrict__ uPtr = addr.upperAddr().begin();
    const label* const __restrict__ lPtr = addr.lowerAddr().begin();
    const label* const __restrict__ ownStartPtr = addr.ownerStartAddr().begin();
    const label* const __restrict__ losortStartAddrPtr =
        addr.losortStartAddr().begin();
    const label* const __restrict__ losortAddrPtr = addr.losortAddr().begin();

    forAll(seedFaces, i)
    {
        allFaceInfo[seedFaces[i]] = seedInfo[i];
    }

    DynamicList<label> changedFaces(seedFaces);
    DynamicList<label> changedCells;
    while (true)
    {
        // face to cell
        changedCells.clear();
        for (const label facei : changedFaces)
        {
            const point& origin = allFaceInfo[facei];
            const point& fc = faceCentres[facei];

            if (closer(fc, origin, allCellInfo[lPtr[facei]]))
            {
                changedCells.append(lPtr[facei]);
            }
            if (closer(fc, origin, allCellInfo[uPtr[facei]]))
            {
                changedCells.append(uPtr[facei]);
            }
        }


        // cell to face
        changedFaces.clear();
        for (const label celli : changedCells)
        {
            const auto& origin = allCellInfo[celli];
            const point& cc = cellCentres[celli];

            {
                const label fStart = ownStartPtr[celli];
                const label fEnd = ownStartPtr[celli + 1];

                for (label facei=fStart; facei<fEnd; facei++)
                {
                    if (closer(cc, origin, allFaceInfo[facei]))
                    {
                        changedFaces.append(facei);
                    }
                }
            }
            {
                const label fStart = losortStartAddrPtr[celli];
                const label fEnd = losortStartAddrPtr[celli + 1];

                for (label i=fStart; i<fEnd; i++)
                {
                    label facei = losortAddrPtr[i];
                    if (closer(cc, origin, allFaceInfo[facei]))
                    {
                        changedFaces.append(facei);
                    }
                }
            }
        }

        if (changedFaces.empty())
        {
            break;
        }
    }
}




void writeAgglomeration
(
    const bool normalise,
    const fvMesh& mesh,
    const labelList& cellToBaseMesh
)
{
    volScalarField scalarAgglomeration
    (
        IOobject
        (
            "agglomeration",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    scalarField& fld = scalarAgglomeration.primitiveFieldRef();
    forAll(fld, celli)
    {
        fld[celli] = cellToBaseMesh[celli];
    }
    if (normalise)
    {
        fld /= max(fld);
    }
    scalarAgglomeration.correctBoundaryConditions();
    scalarAgglomeration.write();
}


// Main program:

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "normalise",
        "normalise agglomeration (0..1)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool normalise = args.found("normalise");

    #include "createMesh.H"

    const fvSolution& sol = static_cast<const fvSolution&>(mesh);
    const dictionary& pDict = sol.subDict("solvers").subDict("p");

    const GAMGAgglomeration& agglom = GAMGAgglomeration::New
    (
        mesh,
        pDict
    );

    // Current level geometry
    scalarField cellVolumes(mesh.cellVolumes());
    pointField cellCentres(mesh.cellCentres());
    scalarField faceAreas(mag(mesh.faceAreas()));
    pointField faceCentres(mesh.faceCentres());

    // Mapping back to base (= fvMesh)
    labelList cellToBase(identity(mesh.nCells()));
    labelListList baseToCell(invertOneToMany(mesh.nCells(), cellToBase));

    // Write initial agglomeration
    Info<< "Writing initial cell distribution to "
        << runTime.timeName() << endl;
    writeAgglomeration(false, mesh, cellToBase);

    for (label level = 0; level < agglom.size(); level++)
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const lduMesh& fineMesh = agglom.meshLevel(level);
        const labelList& addr = agglom.restrictAddressing(level);
        label coarseSize = max(addr)+1;

        Info<< "Level : " << level << endl
            << "    current size      : "
            << returnReduce(addr.size(), sumOp<label>()) << endl
            << "    agglomerated size : "
            << returnReduce(coarseSize, sumOp<label>()) << endl;



        // Update adressing to/from base mesh
        forAll(addr, fineI)
        {
            const labelList& cellLabels = baseToCell[fineI];
            forAll(cellLabels, i)
            {
                cellToBase[cellLabels[i]] = addr[fineI];
            }
        }
        baseToCell = invertOneToMany(coarseSize, cellToBase);


        // Write agglomeration
        Info<< "Writing cell distribution to "
            << runTime.timeName() << endl;
        writeAgglomeration(normalise, mesh, cellToBase);



        // Write lines inbetween cellCentres
        {
            const auto& addr = fineMesh.lduAddr();
            const auto& lower = addr.lowerAddr();
            const auto& upper = addr.upperAddr();

            OBJstream os(runTime.timePath()/"agglomeration.obj");
            forAll(lower, facei)
            {
                const point& lowerCc = cellCentres[lower[facei]];
                const point& upperCc = cellCentres[upper[facei]];
                os.write(linePointRef(lowerCc, upperCc));
            }
        }
        

        // Update volumes, centroids for coarse level
        {
            scalarField coarseVolumes(coarseSize);
            pointField coarseCentroids(coarseSize);

            // Agglomerate/restrict
            coarseVolumes = Zero;
            coarseCentroids = Zero;
            forAll(addr, i)
            {
                coarseVolumes[addr[i]] += cellVolumes[i];
                coarseCentroids[addr[i]] += cellVolumes[i]*cellCentres[i];
            }
            coarseCentroids /= coarseVolumes;

            cellVolumes = std::move(coarseVolumes);
            cellCentres = std::move(coarseCentroids);
        }
        {
            scalarField coarseAreas(coarseSize);
            pointField coarseCentroids(coarseSize);

            // Agglomerate/restrict
            coarseAreas = Zero;
            coarseCentroids = Zero;
            forAll(addr, i)
            {
                coarseAreas[addr[i]] += faceAreas[i];
                coarseCentroids[addr[i]] += faceAreas[i]*faceCentres[i];
            }
            coarseCentroids /= coarseAreas;

            faceAreas = std::move(coarseAreas);
            faceCentres = std::move(coarseCentroids);
        }

        Info<< endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
