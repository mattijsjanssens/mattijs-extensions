/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    snappyHexMesh

Description
    Automatic split hex mesher. Refines and snaps to surface.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "snappyRefineDriver.H"
#include "snappySnapDriver.H"
#include "snappyLayerDriver.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "decompositionMethod.H"
#include "noDecomp.H"
#include "fvMeshDistribute.H"
#include "wallPolyPatch.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "layerParameters.H"
#include "vtkSetWriter.H"
#include "faceSet.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "cellModeller.H"
#include "uindirectPrimitivePatch.H"
#include "surfZoneIdentifierList.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurface.H"
#include "globalIndex.H"
#include "IOmanip.H"
#include "fvMeshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convert size (as fraction of defaultCellSize) to refinement level
label sizeCoeffToRefinement
(
    const scalar level0Coeff,   // ratio of hex cell size v.s. defaultCellSize
    const scalar sizeCoeff
)
{
     return round(::log(level0Coeff/sizeCoeff)/::log(2));
}


autoPtr<refinementSurfaces> createRefinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& surfacesDict,
    const dictionary& shapeControlDict,
    const label gapLevelIncrement,
    const scalar level0Coeff
)
{
    autoPtr<refinementSurfaces> surfacePtr;

    // Count number of surfaces.
    label surfI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry.names()[geomI];

        if (surfacesDict.found(geomName))
        {
            surfI++;
        }
    }

    labelList surfaces(surfI);
    wordList names(surfI);
    PtrList<surfaceZonesInfo> surfZones(surfI);

    labelList regionOffset(surfI);

    labelList globalMinLevel(surfI, 0);
    labelList globalMaxLevel(surfI, 0);
    labelList globalLevelIncr(surfI, 0);
    PtrList<dictionary> globalPatchInfo(surfI);
    List<Map<label>> regionMinLevel(surfI);
    List<Map<label>> regionMaxLevel(surfI);
    List<Map<label>> regionLevelIncr(surfI);
    List<Map<scalar>> regionAngle(surfI);
    List<Map<autoPtr<dictionary>>> regionPatchInfo(surfI);

    HashSet<word> unmatchedKeys(surfacesDict.toc());

    surfI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry.names()[geomI];

        const entry* ePtr = surfacesDict.lookupEntryPtr(geomName, false, true);

        if (ePtr)
        {
            const dictionary& shapeDict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            names[surfI] = geomName;
            surfaces[surfI] = geomI;

            const searchableSurface& surface = allGeometry[geomI];

            // Find the index in shapeControlDict
            // Invert surfaceCellSize to get the refinementLevel

            const word scsFuncName =
                shapeDict.lookup("surfaceCellSizeFunction");
            const dictionary& scsDict =
                shapeDict.subDict(scsFuncName + "Coeffs");

            const scalar surfaceCellSize =
                readScalar(scsDict.lookup("surfaceCellSizeCoeff"));

            const label refLevel = sizeCoeffToRefinement
            (
                level0Coeff,
                surfaceCellSize
            );

            globalMinLevel[surfI] = refLevel;
            globalMaxLevel[surfI] = refLevel;
            globalLevelIncr[surfI] = gapLevelIncrement;

            // Surface zones
            surfZones.set(surfI, new surfaceZonesInfo(surface, shapeDict));


            // Global perpendicular angle
            if (shapeDict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfI,
                    shapeDict.subDict("patchInfo").clone()
                );
            }


            // Per region override of patchInfo

            if (shapeDict.found("regions"))
            {
                const dictionary& regionsDict = shapeDict.subDict("regions");
                const wordList& regionNames =
                    allGeometry[surfaces[surfI]].regions();

                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfI].insert
                            (
                                regionI,
                                regionDict.subDict("patchInfo").clone()
                            );
                        }
                    }
                }
            }

            // Per region override of cellSize
            if (shapeDict.found("regions"))
            {
                const dictionary& shapeControlRegionsDict =
                    shapeDict.subDict("regions");
                const wordList& regionNames =
                    allGeometry[surfaces[surfI]].regions();

                forAll(regionNames, regionI)
                {
                    if (shapeControlRegionsDict.found(regionNames[regionI]))
                    {
                        const dictionary& shapeControlRegionDict =
                            shapeControlRegionsDict.subDict
                            (
                                regionNames[regionI]
                            );

                        const word scsFuncName =
                            shapeControlRegionDict.lookup
                            (
                                "surfaceCellSizeFunction"
                            );
                        const dictionary& scsDict =
                            shapeControlRegionDict.subDict
                            (
                                scsFuncName + "Coeffs"
                            );

                        const scalar surfaceCellSize =
                            readScalar
                            (
                                scsDict.lookup("surfaceCellSizeCoeff")
                            );

                        const label refLevel = sizeCoeffToRefinement
                        (
                            level0Coeff,
                            surfaceCellSize
                        );

                        regionMinLevel[surfI].insert(regionI, refLevel);
                        regionMaxLevel[surfI].insert(regionI, refLevel);
                        regionLevelIncr[surfI].insert(regionI, 0);
                    }
                }
            }

            surfI++;
        }
    }

    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces, surfI)
    {
        regionOffset[surfI] = nRegions;
        nRegions += allGeometry[surfaces[surfI]].regions().size();
    }

    // Rework surface specific information into information per global region
    labelList minLevel(nRegions, 0);
    labelList maxLevel(nRegions, 0);
    labelList gapLevel(nRegions, -1);
    PtrList<dictionary> patchInfo(nRegions);

    forAll(globalMinLevel, surfI)
    {
        label nRegions = allGeometry[surfaces[surfI]].regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegionI = regionOffset[surfI] + i;
            minLevel[globalRegionI] = globalMinLevel[surfI];
            maxLevel[globalRegionI] = globalMaxLevel[surfI];
            gapLevel[globalRegionI] =
                maxLevel[globalRegionI]
              + globalLevelIncr[surfI];

            if (globalPatchInfo.set(surfI))
            {
                patchInfo.set
                (
                    globalRegionI,
                    globalPatchInfo[surfI].clone()
                );
            }
        }

        // Overwrite with region specific information
        forAllConstIter(Map<label>, regionMinLevel[surfI], iter)
        {
            label globalRegionI = regionOffset[surfI] + iter.key();

            minLevel[globalRegionI] = iter();
            maxLevel[globalRegionI] = regionMaxLevel[surfI][iter.key()];
            gapLevel[globalRegionI] =
                maxLevel[globalRegionI]
              + regionLevelIncr[surfI][iter.key()];
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfI];
        forAllConstIter(Map<autoPtr<dictionary>>, localInfo, iter)
        {
            label globalRegionI = regionOffset[surfI] + iter.key();
            patchInfo.set(globalRegionI, iter()().clone());
        }
    }

    surfacePtr.set
    (
        new refinementSurfaces
        (
            allGeometry,
            surfaces,
            names,
            surfZones,
            regionOffset,
            minLevel,
            maxLevel,
            gapLevel,
            scalarField(nRegions, -GREAT),  //perpendicularAngle,
            patchInfo
        )
    );


    const refinementSurfaces& rf = surfacePtr();

    // Determine maximum region name length
    label maxLen = 0;
    forAll(rf.surfaces(), surfI)
    {
        label geomI = rf.surfaces()[surfI];
        const wordList& regionNames = allGeometry.regionNames()[geomI];
        forAll(regionNames, regionI)
        {
            maxLen = Foam::max(maxLen, label(regionNames[regionI].size()));
        }
    }


    Info<< setw(maxLen) << "Region"
        << setw(10) << "Min Level"
        << setw(10) << "Max Level"
        << setw(10) << "Gap Level" << nl
        << setw(maxLen) << "------"
        << setw(10) << "---------"
        << setw(10) << "---------"
        << setw(10) << "---------" << endl;

    forAll(rf.surfaces(), surfI)
    {
        label geomI = rf.surfaces()[surfI];

        Info<< rf.names()[surfI] << ':' << nl;

        const wordList& regionNames = allGeometry.regionNames()[geomI];

        forAll(regionNames, regionI)
        {
            label globalI = rf.globalRegion(surfI, regionI);

            Info<< setw(maxLen) << regionNames[regionI]
                << setw(10) << rf.minLevel()[globalI]
                << setw(10) << rf.maxLevel()[globalI]
                << setw(10) << rf.gapLevel()[globalI] << endl;
        }
    }


    return surfacePtr;
}


void extractSurface
(
    const polyMesh& mesh,
    const Time& runTime,
    const labelHashSet& includePatches,
    const fileName& outFileName
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // Collect sizes. Hash on names to handle local-only patches (e.g.
    //  processor patches)
    HashTable<label> patchSize(1000);
    label nFaces = 0;
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        patchSize.insert(pp.name(), pp.size());
        nFaces += pp.size();
    }
    Pstream::mapCombineGather(patchSize, plusEqOp<label>());


    // Allocate zone/patch for all patches
    HashTable<label> compactZoneID(1000);
    forAllConstIter(HashTable<label>, patchSize, iter)
    {
        label sz = compactZoneID.size();
        compactZoneID.insert(iter.key(), sz);
    }
    Pstream::mapCombineScatter(compactZoneID);


    // Rework HashTable into labelList just for speed of conversion
    labelList patchToCompactZone(bMesh.size(), -1);
    forAllConstIter(HashTable<label>, compactZoneID, iter)
    {
        label patchi = bMesh.findPatchID(iter.key());
        if (patchi != -1)
        {
            patchToCompactZone[patchi] = iter();
        }
    }

    // Collect faces on zones
    DynamicList<label> faceLabels(nFaces);
    DynamicList<label> compactZones(nFaces);
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        forAll(pp, i)
        {
            faceLabels.append(pp.start()+i);
            compactZones.append(patchToCompactZone[pp.index()]);
        }
    }

    // Addressing engine for all faces
    uindirectPrimitivePatch allBoundary
    (
        UIndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );


    // Find correspondence to master points
    labelList pointToGlobal;
    labelList uniqueMeshPoints;
    autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
    (
        allBoundary.meshPoints(),
        allBoundary.meshPointMap(),
        pointToGlobal,
        uniqueMeshPoints
    );

    // Gather all unique points on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()] = pointField
    (
        mesh.points(),
        uniqueMeshPoints
    );
    Pstream::gatherList(gatheredPoints);

    // Gather all faces
    List<faceList> gatheredFaces(Pstream::nProcs());
    gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
    forAll(gatheredFaces[Pstream::myProcNo()], i)
    {
        inplaceRenumber(pointToGlobal, gatheredFaces[Pstream::myProcNo()][i]);
    }
    Pstream::gatherList(gatheredFaces);

    // Gather all ZoneIDs
    List<labelList> gatheredZones(Pstream::nProcs());
    gatheredZones[Pstream::myProcNo()] = compactZones.xfer();
    Pstream::gatherList(gatheredZones);

    // On master combine all points, faces, zones
    if (Pstream::master())
    {
        pointField allPoints = ListListOps::combine<pointField>
        (
            gatheredPoints,
            accessOp<pointField>()
        );
        gatheredPoints.clear();

        faceList allFaces = ListListOps::combine<faceList>
        (
            gatheredFaces,
            accessOp<faceList>()
        );
        gatheredFaces.clear();

        labelList allZones = ListListOps::combine<labelList>
        (
            gatheredZones,
            accessOp<labelList>()
        );
        gatheredZones.clear();


        // Zones
        surfZoneIdentifierList surfZones(compactZoneID.size());
        forAllConstIter(HashTable<label>, compactZoneID, iter)
        {
            surfZones[iter()] = surfZoneIdentifier(iter.key(), iter());
            Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
                << endl;
        }

        UnsortedMeshedSurface<face> unsortedFace
        (
            xferMove(allPoints),
            xferMove(allFaces),
            xferMove(allZones),
            xferMove(surfZones)
        );


        MeshedSurface<face> sortedFace(unsortedFace);

        fileName globalCasePath
        (
            runTime.processorCase()
          ? runTime.path()/".."/outFileName
          : runTime.path()/outFileName
        );

        Info<< "Writing merged surface to " << globalCasePath << endl;

        sortedFace.write(globalCasePath);
    }
}


// Check writing tolerance before doing any serious work
scalar getMergeDistance(const polyMesh& mesh, const scalar mergeTol)
{
    const boundBox& meshBb = mesh.bounds();
    scalar mergeDist = mergeTol * meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    // check writing tolerance
    if (mesh.time().writeFormat() == IOstream::ASCII)
    {
        const scalar writeTol = std::pow
        (
            scalar(10.0),
            -scalar(IOstream::defaultPrecision())
        );

        if (mergeTol < writeTol)
        {
            FatalErrorInFunction
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << nl
                << "Your merging tolerance (" << mergeTol
                << ") is finer than this." << nl
                << "Change to binary writeFormat, "
                << "or increase the writePrecision" << endl
                << "or adjust the merge tolerance (mergeTol)."
                << exit(FatalError);
        }
    }

    return mergeDist;
}


void removeZeroSizedPatches(fvMesh& mesh)
{
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduce(pp.size(), sumOp<label>())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
        Info<< endl;
    }
}


// Write mesh and additional information
void writeMesh
(
    const string& msg,
    const meshRefinement& meshRefiner,
    const meshRefinement::debugType debugLevel,
    const meshRefinement::writeType writeLevel
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    meshRefiner.printMeshInfo(debugLevel, msg);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    //label flag = meshRefinement::MESH;
    //if (writeLevel)
    //{
    //    flag |= meshRefinement::SCALARLEVELS;
    //}
    //if (debug & meshRefinement::OBJINTERSECTIONS)
    //{
    //    flag |= meshRefinement::OBJINTERSECTIONS;
    //}
    meshRefiner.write
    (
        debugLevel,
        meshRefinement::writeType(writeLevel | meshRefinement::WRITEMESH),
        mesh.time().path()/meshRefiner.timeName()
    );
    Info<< "Wrote mesh in = "
        << mesh.time().cpuTimeIncrement() << " s." << endl;
}


int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    // Check patches and faceZones are synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    meshRefinement::checkCoupledFaceZones(mesh);


    // Read meshing dictionary
    const word dictName("snappyHexMeshDict");
    #include "setSystemMeshDictionaryIO.H"
    const IOdictionary meshDict(dictIO);

    // all surface geometry
    const dictionary& geometryDict = meshDict.subDict("geometry");

    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            //mesh.time().findInstance("triSurface", word::null),// instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        true    //meshDict.lookupOrDefault("singleRegionName", true)
    );


    forAll(allGeometry, geomi)
    {
        Pout<< "geom:" << allGeometry[geomi].name()
            << " size:" << allGeometry[geomi].size()
            << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
