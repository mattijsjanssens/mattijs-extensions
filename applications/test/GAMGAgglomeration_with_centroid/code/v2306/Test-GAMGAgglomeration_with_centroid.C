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
        if (fineVals[i] != nullValue)
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
        if (fineVals[i] != nullValue)
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
        if (fineVals[i] != nullValue)
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

    // Current level geometry
    scalarField cellVolumes(mesh.cellVolumes());
    pointField cellCentres(mesh.cellCentres());
    scalarField faceAreas(mag(mesh.faceAreas()));
    pointField faceCentres(mesh.faceCentres());

    // Nearest wall info
    scalarField cellWallDistance(mesh.nCells(), GREAT);
    pointField cellNearestWall(mesh.nCells(), point::uniform(GREAT));
    for (const label patchi : patchIDs)
    {
        const auto& fvp = mesh.boundary()[patchi];
        const auto& faceCells = fvp.faceCells();
        UIndirectList<point>(cellNearestWall, faceCells) =
            mesh.C().boundaryField()[patchi];

        forAll(faceCells, i)
        {
            const label celli = faceCells[i];
            cellWallDistance[celli] =
                Foam::mag(cellCentres[celli]-cellNearestWall[celli]);
        }
    }


    // Per base cell the agglomerated cell
    labelList baseToCell(identity(mesh.nCells()));

    // Write initial agglomeration
    Info<< "Writing initial cell distribution to "
        << runTime.timeName() << endl;
    writeAgglomeration(false, mesh, baseToCell);

    for (label level = 0; level < agglom.size(); level++)
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const lduMesh& fineMesh = agglom.meshLevel(level);
        const auto& fineAddr = fineMesh.lduAddr();
        const label nFineCells = fineAddr.size();
        const label nFineFaces = fineAddr.lowerAddr().size();
        const labelList& fineToCoarse = agglom.restrictAddressing(level);
        const label coarseSize = max(fineToCoarse)+1;

        Info<< "Level : " << level << endl
            << "    current size      : "
            << returnReduce(fineToCoarse.size(), sumOp<label>()) << endl
            << "    agglomerated size : "
            << returnReduce(coarseSize, sumOp<label>()) << endl;

        if (max(baseToCell)+1 != fineToCoarse.size())
        {
            FatalErrorInFunction << "baseToCell:" << max(baseToCell)+1
                << " fineToCoarse.size():" << fineToCoarse.size()
                << exit(FatalError);
        }
        if
        (
            fineToCoarse.size() != nFineCells
         || cellCentres.size() != nFineCells
        )
        {
            FatalErrorInFunction << exit(FatalError);
        }

        //Override faceCentres since agglomerated face centres not
        // representative
        faceCentres = interpolate(fineMesh, cellCentres);


        // For postprocessing: mapped near-wall distance
        Info<< "Writing cellWallDistance to " << mesh.time().timeName()
            << endl;
        writeField("cellWallDistance", mesh, baseToCell, cellWallDistance);


        {
            List<point> allCellInfo(fineToCoarse.size(), point::uniform(GREAT));
            List<point> allFaceInfo(nFineFaces, point::uniform(GREAT));

            OBJstream os(runTime.timePath()/"seedCells.obj");
            DynamicList<label> seedCells;
            DynamicList<point> seedInfo;
            forAll(cellNearestWall, celli)
            {
                if (cellNearestWall[celli] != point::uniform(GREAT))
                {
                    seedCells.append(celli);
                    seedInfo.append(cellNearestWall[celli]);

                    const point& cc = cellCentres[celli];
                    const point& wallPt = cellNearestWall[celli];
                    os.write(linePointRef(cc, wallPt));
                }
            }
            Pout<< "** seeding " << seedCells.size()
                << " out of " << cellNearestWall.size()
                << endl;

            meshWaveNearest
            (
                allCellInfo.size()+1,   // maxIter
                fineMesh,
                faceCentres,
                cellCentres,
                seedCells,
                seedInfo,
                allFaceInfo,
                allCellInfo
            );


            writeLines
            (
                runTime.timePath()/"allCellInfo.obj",
                cellCentres,
                allCellInfo
            );

            Info<< "Writing allCellInfo to " << mesh.time().timeName() << endl;
            writeField
            (
                "allCellInfo",
                mesh,
                baseToCell,
                mag(cellCentres-allCellInfo)
            );
        }

        // Write agglomeration
        Info<< "Writing cell distribution to "
            << runTime.timeName() << endl;
        writeAgglomeration(normalise, mesh, baseToCell);


        // Write lines inbetween cellCentres
        writeLines
        (
            runTime.timePath()/"agglomeration.obj",
            pointList(cellCentres, fineAddr.lowerAddr()),
            pointList(cellCentres, fineAddr.upperAddr())
        );


        // Calculate normalised volume weights
        scalarField volWeights;
        {
            scalarField sumWeights(cellVolumes.size(), Zero);
            forAll(fineToCoarse, i)
            {
                sumWeights[fineToCoarse[i]] += cellVolumes[i];
            }
            volWeights = cellVolumes;
            forAll(fineToCoarse, i)
            {
                volWeights[i] /= sumWeights[fineToCoarse[i]];
            }
        }
        // Calculate normalised face-area weights
        scalarField faceWeights;
        {
            scalarField sumWeights(faceAreas.size(), Zero);
            forAll(fineToCoarse, i)
            {
                sumWeights[fineToCoarse[i]] += faceAreas[i];
            }
            faceWeights = faceAreas;
            forAll(fineToCoarse, i)
            {
                faceWeights[i] /= sumWeights[fineToCoarse[i]];
            }
        }


        // Update volumes, centroids for coarse level
        cellVolumes = restrictField
        (
            coarseSize,
            fineToCoarse,
            cellVolumes,
            GREAT
        );
        cellCentres = restrictField
        (
            coarseSize,
            fineToCoarse,
            volWeights,
            cellCentres,
            point::uniform(GREAT)
        );

        // Update face geometry for coarse level
        faceAreas = restrictField(coarseSize, fineToCoarse, faceAreas, GREAT);
        faceCentres = restrictField
        (
            coarseSize,
            fineToCoarse,
            faceWeights,
            faceCentres,
            point::uniform(GREAT)
        );

        // Update nearest
        cellWallDistance = restrictUnnormalisedField    //restrictField
        (
            coarseSize,
            fineToCoarse,
            volWeights,
            cellWallDistance,
            GREAT
        );
        cellNearestWall = restrictUnnormalisedField     //restrictField
        (
            coarseSize,
            fineToCoarse,
            volWeights,
            cellNearestWall,
            point::uniform(GREAT)
        );

        // Update adressing to/from base mesh
        forAll(baseToCell, baseCelli)
        {
            baseToCell[baseCelli] = fineToCoarse[baseToCell[baseCelli]];
        }

        Info<< endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
