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
#include "wallDistAddressing.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>> restrictField
(
    const label nCoarse,
    const labelUList& fineToCoarse,
    const Field<Type>& fineVals,
    const Type& nullValue
)
{
    tmp<Field<Type>> tfld(tmp<Field<Type>>::New(nCoarse, Zero));
    Field<Type>& fld = tfld.ref();

    // Agglomerate/restrict
    forAll(fineToCoarse, i)
    {
        if (fineToCoarse[i] >= 0 && fineVals[i] != nullValue)
        {
            fld[fineToCoarse[i]] += fineVals[i];
        }
    }
    return tfld;
}


template<class Type>
tmp<Field<Type>> restrictField
(
    const label nCoarse,
    const labelUList& fineToCoarse,
    const scalarField& fineWeights,
    const Field<Type>& fineVals,
    const Type& nullValue
)
{
    tmp<Field<Type>> tfld(tmp<Field<Type>>::New(nCoarse, Zero));
    Field<Type>& fld = tfld.ref();

    // Agglomerate/restrict
    forAll(fineToCoarse, i)
    {
        if (fineToCoarse[i] >= 0 && fineVals[i] != nullValue)
        {
            fld[fineToCoarse[i]] += fineWeights[i]*fineVals[i];
        }
    }
    return tfld;
}


template<class Type>
tmp<Field<Type>> restrictUnnormalisedField
(
    const label nCoarse,
    const labelUList& fineToCoarse,
    const scalarField& fineWeights,
    const Field<Type>& fineVals,
    const Type& nullValue
)
{
    tmp<Field<Type>> tfld(tmp<Field<Type>>::New(nCoarse, Zero));
    Field<Type>& fld = tfld.ref();
    scalarField coarseSumWeights(nCoarse, Zero);

    // Agglomerate/restrict
    forAll(fineToCoarse, i)
    {
        if (fineToCoarse[i] >= 0 && fineVals[i] != nullValue)
        {
            fld[fineToCoarse[i]] += fineWeights[i]*fineVals[i];
            coarseSumWeights[fineToCoarse[i]] += fineWeights[i];
        }
    }

    forAll(fld, coarsei)
    {
        if (coarseSumWeights[coarsei] == scalar(Zero))
        {
            fld[coarsei] = nullValue;
        }
        else
        {
            fld[coarsei] /= coarseSumWeights[coarsei];
        }
    }

    return tfld;
}


template<class Type>
tmp<Field<Type>> prolongField
(
    const labelUList& fineToCoarse,
    const Field<Type>& coarseVals,
    const Type& nullValue
)
{
    tmp<Field<Type>> tfld(tmp<Field<Type>>::New(fineToCoarse.size(), Zero));
    Field<Type>& fld = tfld.ref();

    // Agglomerate/restrict
    forAll(fineToCoarse, i)
    {
        const auto& coarseVal = coarseVals[fineToCoarse[i]];
        if (coarseVal != nullValue)
        {
            fld[i] = coarseVal;
        }
    }
    return tfld;
}
template<class Type>
void prolongAndAdd
(
    const labelUList& fineToCoarse,
    const Field<Type>& coarseVals,
    const Type& nullValue,
    List<Type>& result
)
{
    tmp<Field<Type>> tfld
    (
        prolongField
        (
            fineToCoarse,
            coarseVals,
            point::uniform(GREAT)
        )
    );
    const auto& fld = tfld();
    forAll(fld, i)
    {
        if (fld[i] != nullValue)
        {
            result[i] += fld[i];
        }
    }
}


template<class Type>
tmp<Field<Type>> interpolate
(
    const lduMesh& mesh,
    const Field<Type>& vals
)
{
    const lduAddressing& addr = mesh.lduAddr();
    const label nFaces = addr.upperAddr().size();
    const label* const __restrict__ uPtr = addr.upperAddr().begin();
    const label* const __restrict__ lPtr = addr.lowerAddr().begin();

    tmp<Field<Type>> tfld(tmp<Field<Type>>::New(nFaces));
    auto& fld = tfld.ref();

    for (label facei = 0; facei < nFaces; facei++)
    {
        fld[facei] = 0.5*(vals[lPtr[facei]]+vals[uPtr[facei]]);
    }
    return tfld;
}


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


label meshWaveNearest
(
    const label maxIter,
    const lduMesh& mesh,
    const pointField& faceCentres,
    const pointField& cellCentres,
    const labelUList& seedCells,
    const List<point>& seedInfo,
    List<point>& allFaceInfo,
    List<point>& allCellInfo
)
{
    // Walk from cell to face since boundary conditions given as initial
    // cells.


    const lduAddressing& addr = mesh.lduAddr();

    if
    (
        (allCellInfo.size() != addr.size())
     || (allFaceInfo.size() != addr.upperAddr().size())
    )
    {
        FatalErrorInFunction<< "illegal container sizes" << exit(FatalError);
    }

    const label* const __restrict__ uPtr = addr.upperAddr().begin();
    const label* const __restrict__ lPtr = addr.lowerAddr().begin();
    const label* const __restrict__ ownStartPtr = addr.ownerStartAddr().begin();
    const label* const __restrict__ losortStartAddrPtr =
        addr.losortStartAddr().begin();
    const label* const __restrict__ losortAddrPtr = addr.losortAddr().begin();

    Pout<< "starting with seedCells:" << seedCells.size() << endl;

    forAll(seedCells, i)
    {
        allCellInfo[seedCells[i]] = seedInfo[i];
    }


    DynamicList<label> changedCells(seedCells);
    bitSet isChangedCell(addr.size());
    DynamicList<label> changedFaces;
    bitSet isChangedFace(addr.lowerAddr().size());

    label iter = 0;
    for (;iter < maxIter; iter++)
    {
        // cell to face
        changedFaces.clear();
        isChangedFace = false;
        
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
                        if (isChangedFace.set(facei))
                        {
                            changedFaces.append(facei);
                        }
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
                        if (isChangedFace.set(facei))
                        {
                            changedFaces.append(facei);
                        }
                    }
                }
            }
        }
        Pout<< "   from changedCells:" << changedCells.size()
            << " to changedFaces:" << changedFaces.size() << endl;

        // face to cell
        changedCells.clear();
        isChangedCell = false;
        for (const label facei : changedFaces)
        {
            const point& origin = allFaceInfo[facei];
            const point& fc = faceCentres[facei];

            if (closer(fc, origin, allCellInfo[lPtr[facei]]))
            {
                if (isChangedCell.set(lPtr[facei]))
                {
                    changedCells.append(lPtr[facei]);
                }
            }
            if (closer(fc, origin, allCellInfo[uPtr[facei]]))
            {
                if (isChangedCell.set(uPtr[facei]))
                {
                    changedCells.append(uPtr[facei]);
                }
            }
        }
        Pout<< "   from changedFaces:" << changedFaces.size()
            << " to changedCells:" << changedCells.size() << endl;

        if (changedCells.empty())
        {
            break;
        }
    }
    return iter;
}
void solveWallDistance
(
    const label maxIter,
    const lduMesh& mesh,
    const pointField& cellNearestWall,
    const pointField& faceCentres,
    const pointField& cellCentres,
    pointList& allFaceInfo,
    pointList& allCellInfo,

    pointField& allFaceChange,
    pointField& allCellChange
)
{
    // Solve mesh starting from boundaries
    DynamicList<label> seedCells;
    DynamicList<point> seedInfo;
    forAll(cellNearestWall, celli)
    {
        if (cellNearestWall[celli] != point::uniform(GREAT))
        {
            seedCells.append(celli);
            seedInfo.append(cellNearestWall[celli]);
        }
    }
    Pout<< "** seeding " << seedCells.size()
        << " out of " << cellNearestWall.size()
        << endl;


    allFaceChange = allFaceInfo;
    allCellChange = allCellInfo;

    meshWaveNearest
    (
        maxIter,
        mesh,
        faceCentres,
        cellCentres,
        seedCells,
        seedInfo,
        allFaceInfo,
        allCellInfo
    );

    forAll(allFaceInfo, i)
    {
        if (allFaceInfo[i] != point::uniform(GREAT))
        {
            allFaceChange[i] -= allFaceInfo[i];
        }
        else
        {
            allFaceChange[i] = Zero;
        }
    }

    forAll(allCellInfo, i)
    {
        if (allCellInfo[i] != point::uniform(GREAT))
        {
            allCellChange[i] -= allCellInfo[i];
        }
        else
        {
            allCellChange[i] = Zero;
        }
    }
}


void writeField
(
    const word& name,
    const fvMesh& baseMesh,
    const labelList& baseToCell,
    const scalarField& coarseFld
)
{
    volScalarField fld
    (
        IOobject
        (
            name,
            baseMesh.time().timeName(),
            baseMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        baseMesh,
        dimensionedScalar(dimless, Zero)
    );
    forAll(baseToCell, baseCelli)
    {
        const label celli = baseToCell[baseCelli];
        fld[baseCelli] = coarseFld[celli];
    }
    fld.correctBoundaryConditions();
    fld.write();
}


void writeAgglomeration
(
    const bool normalise,
    const fvMesh& mesh,
    const labelList& meshToAgglom
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
        fld[celli] = meshToAgglom[celli];
    }
    if (normalise)
    {
        fld /= max(fld);
    }
    scalarAgglomeration.correctBoundaryConditions();
    scalarAgglomeration.write();
}
void writeLines
(
    const fileName& fName,
    const UList<point>& start,
    const UList<point>& end
)
{
    OBJstream os(fName);
    forAll(start, i)
    {
        os.write(linePointRef(start[i], end[i]));
    }
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

    // Get patchIDs
    const labelList patchIDs
    (
        mesh.boundaryMesh().findPatchIDs<wallPolyPatch>().sortedToc()
    );
    Info<< "** Walls:" << patchIDs << endl;

    const GAMGAgglomeration& agglom = GAMGAgglomeration::New
    (
        mesh,
        pDict
    );

    const label nLevels = agglom.size();


    // Current level geometry
    List<scalarField> cellVolumes(nLevels);
    List<pointField> cellCentres(nLevels);
    List<scalarField> faceAreas(nLevels);
    List<pointField> faceCentres(nLevels);
    // Per base cell the agglomerated cell (for postprocessing)
    labelListList baseToCell(nLevels);
    // Nearest wall info
    List<scalarField> cellWallDistance(nLevels);
    List<pointField> cellNearestWall(nLevels);
    // Current solved-for variables
    List<pointList> allCellInfo(nLevels);
    List<pointField> allCellResidual(nLevels);
    List<pointList> allFaceInfo(nLevels);
    List<pointField> allFaceResidual(nLevels);


    // Fill level 0 from the current mesh
    {
        cellVolumes[0] = mesh.cellVolumes();
        cellCentres[0] = mesh.cellCentres();
        faceAreas[0] = mag(mesh.faceAreas());
        faceCentres[0] = mesh.faceCentres();

        scalarField& wallDist = cellWallDistance[0];
        pointField& nearestWall = cellNearestWall[0];
        wallDist.setSize(mesh.nCells(), GREAT);
        nearestWall.setSize(mesh.nCells(), point::uniform(GREAT));
        for (const label patchi : patchIDs)
        {
            const auto& fvp = mesh.boundary()[patchi];
            const auto& faceCells = fvp.faceCells();
            const auto& fC = mesh.C().boundaryField()[patchi];

            forAll(faceCells, i)
            {
                const label celli = faceCells[i];
                nearestWall[celli] = fC[i];
                wallDist[celli] = Foam::mag(mesh.C()[celli]-fC[i]);
            }
        }
        baseToCell[0] = identity(mesh.nCells());


        allCellInfo[0] = pointList(mesh.nCells(), point::uniform(GREAT));
        allFaceInfo[0] = pointList(mesh.nFaces(), point::uniform(GREAT));

        // Write initial agglomeration
        Info<< "Writing initial cell distribution to "
            << runTime.timeName() << endl;
        writeAgglomeration(false, mesh, baseToCell[0]);
    }

    // Fill other levels
    for (label leveli = 1; leveli < nLevels; leveli++)
    {
        const lduMesh& levelMesh = agglom.meshLevel(leveli);
        const auto& levelAddr = levelMesh.lduAddr();

        // Addressing from previous to current mesh
        const labelList& fineToCoarseCells =
            agglom.restrictAddressing(leveli-1);
        const labelList& fineToCoarseFaces =
            agglom.faceRestrictAddressing(leveli-1);

        //Pout<< "nFineCells:" << levelAddr.size() << endl;
        //Pout<< "fineToCoarseCells:" << flatOutput(fineToCoarseCells) << endl;
        //const label nCoarseCells = agglom.nCells(leveli);
        //Pout<< "nCoarseCells:" << nCoarseCells << endl;

        //Pout<< "nFineFaces:" << levelAddr.lowerAddr().size() << endl;
        //Pout<< "fineToCoarseFaces:" << flatOutput(fineToCoarseFaces) << endl;
        //const label nCoarseFaces = agglom.nFaces(leveli);
        //Pout<< "nCoarseFaces:" << nCoarseFaces << endl;


        // Update volumes, centroids for coarse level
        cellVolumes[leveli] = restrictField
        (
            levelAddr.size(),
            fineToCoarseCells,
            cellVolumes[leveli-1],
            GREAT
        );

        // Calculate normalised volume weights
        scalarField volWeights(fineToCoarseCells.size());
        {
            forAll(fineToCoarseCells, i)
            {
                volWeights[i] =
                    cellVolumes[leveli-1][i]
                  / cellVolumes[leveli][fineToCoarseCells[i]];
            }
        }

        // Update face geometry for coarse level
        faceAreas[leveli] = restrictField
        (
            levelAddr.lowerAddr().size(),
            fineToCoarseFaces,
            faceAreas[leveli-1],
            GREAT
        );

        // Calculate normalised face-area weights
        scalarField faceWeights(fineToCoarseFaces.size());
        {
            forAll(fineToCoarseFaces, i)
            {
                const label coarsei = fineToCoarseFaces[i];

                if (coarsei >= 0)
                {
                    faceWeights[i] =
                        faceAreas[leveli-1][i]
                      / faceAreas[leveli][coarsei];
                }
            }
        }

        cellCentres[leveli] = restrictField
        (
            levelAddr.size(),
            fineToCoarseCells,
            volWeights,
            cellCentres[leveli-1],
            point::uniform(GREAT)
        );

        //faceCentres[leveli] = restrictField
        //(
        //    levelAddr.lowerAddr().size(),
        //    fineToCoarseFaces,
        //    faceWeights,
        //    faceCentres[leveli-1],
        //    point::uniform(GREAT)
        //);
        // Override faceCentres since agglomerated face centres not
        // representative
        faceCentres[leveli] = interpolate(levelMesh, cellCentres[leveli]);


        // Update nearest
        cellWallDistance[leveli] = restrictUnnormalisedField    //restrictField
        (
            levelAddr.size(),
            fineToCoarseCells,
            volWeights,
            cellWallDistance[leveli-1],
            GREAT
        );
        cellNearestWall[leveli] = restrictUnnormalisedField     //restrictField
        (
            levelAddr.size(),
            fineToCoarseCells,
            volWeights,
            cellNearestWall[leveli-1],
            point::uniform(GREAT)
        );

        // Update adressing to/from base mesh
        baseToCell[leveli] = UIndirectList<label>
        (
            fineToCoarseCells,
            baseToCell[leveli-1]
        );

        allCellInfo[leveli] =
            pointList(levelAddr.size(), point::uniform(GREAT));
        allFaceInfo[leveli] =
            pointList(levelAddr.lowerAddr().size(), point::uniform(GREAT));

        // Write agglomeration
        Info<< "Writing cell distribution to "
            << runTime.timeName() << endl;
        writeAgglomeration(normalise, mesh, baseToCell[leveli]);
    }


    // Write fields
    forAll(baseToCell, leveli)
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const auto& map = baseToCell[leveli];

        // Current level geometry
        writeField("cellVolumes", mesh, map, cellVolumes[leveli]);
        writeField("faceAreas", mesh, map, faceAreas[leveli]);
        // Nearest wall info
        const auto& wallDist = cellWallDistance[leveli];
        writeField("cellWallDistance", mesh, map, wallDist);
        //const auto& nearestWall = cellNearestWall[leveli];
        //writeField("cellNearestWall", mesh, map, nearestWall);
    }


    // V cycles
    for (label iter = 0; iter < 1; iter++)
    {
        // Determine wall distance on the coarsest mesh accurately
        {
            const label leveli = nLevels-1;
            const lduMesh& levelMesh = agglom.meshLevel(leveli);
            const auto& levelAddr = levelMesh.lduAddr();

            solveWallDistance
            (
                levelAddr.size()+1, // maxIter
                levelMesh,
                cellNearestWall[leveli],
                faceCentres[leveli],
                cellCentres[leveli],
                allFaceInfo[leveli],
                allCellInfo[leveli],

                allFaceResidual[leveli],
                allCellResidual[leveli]
            );
        }


        for (label leveli=nLevels-2; leveli>=0; leveli--)
        {
            const lduMesh& levelMesh = agglom.meshLevel(leveli);

            // Map coarse level residual/changes to this level
            prolongAndAdd
            (
                agglom.faceRestrictAddressing(leveli),
                allFaceResidual[leveli+1],
                point::uniform(GREAT),
                allFaceInfo[leveli]
            );
            prolongAndAdd
            (
                agglom.restrictAddressing(leveli),
                allCellResidual[leveli+1],
                point::uniform(GREAT),
                allCellInfo[leveli]
            );

            // Solve a bit less
            solveWallDistance
            (
                10,                 // maxIter
                levelMesh,
                cellNearestWall[leveli],
                faceCentres[leveli],
                cellCentres[leveli],
                allFaceInfo[leveli],
                allCellInfo[leveli],

                allFaceResidual[leveli],
                allCellResidual[leveli]
            );
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
