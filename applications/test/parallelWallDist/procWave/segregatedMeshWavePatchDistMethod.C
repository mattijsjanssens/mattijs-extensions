/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "segregatedMeshWavePatchDistMethod.H"
#include "fvMesh.H"
#include "volFields.H"
//#include "patchDataWave.H"
#include "FaceCellWave.H"
#include "wallPoint.H"
#include "emptyFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(segregatedMeshWave, 0);
    addToRunTimeSelectionTable(patchDistMethod, segregatedMeshWave, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::segregatedMeshWave::segregatedMeshWave
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs),
    correctWalls_(dict.lookupOrDefault<Switch>("correctWalls", true)),
    nUnset_(0)
{}


Foam::patchDistMethods::segregatedMeshWave::segregatedMeshWave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const bool correctWalls
)
:
    patchDistMethod(mesh, patchIDs),
    correctWalls_(correctWalls),
    nUnset_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::segregatedMeshWave::correct(volScalarField& y)
{
    y = dimensionedScalar("yWall", dimLength, great);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();


    //- Wall information for all faces
    List<wallPoint> allFaceInfo(mesh_.nFaces());

    //- Wall information for all cells
    List<wallPoint> allCellInfo(mesh_.nCells());

    //- Pass1: do processor-local distance to wall
Info<< "** PASS1 " << endl;
    {
        label nPatch = 0;
        forAll(pbm, patchi)
        {
            if (patchIDs_.found(patchi))
            {
                nPatch += pbm[patchi].size();
            }
        }

        // Collect patch faces and (one of) its points

        List<wallPoint> faceDist(nPatch);
        labelList changedFaces(nPatch);

        nPatch = 0;
        forAll(pbm, patchi)
        {
            if (patchIDs_.found(patchi))
            {
                const polyPatch& patch = pbm[patchi];
                const vectorField::subField fc = patch.faceCentres();

                // Collect all patch faces
                forAll(fc, patchFacei)
                {
                    label meshFacei = patch.start() + patchFacei;
                    changedFaces[nPatch] = meshFacei;
                    faceDist[nPatch] = wallPoint(fc[patchFacei], 0.0);

                    nPatch++;
                }
            }
        }

        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;

        //- Wave calculation engine.
        FaceCellWave<wallPoint> calc
        (
            mesh_,
            changedFaces,
            faceDist,
            allFaceInfo,
            allCellInfo,
            mesh_.nCells()+1    //maxIter,
        );
        Pstream::parRun() = oldParRun; 
    }


    // Pass2: seed all the processor faces
Info<< "** PASS2 " << endl;
    {
        label nPatch = 0;
        forAll(pbm, patchi)
        {
            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                nPatch += pbm[patchi].size();
            }
        }

        List<wallPoint> faceDist(nPatch);
        labelList changedFaces(nPatch);

        nPatch = 0;
        forAll(pbm, patchi)
        {
            if (isA<processorPolyPatch>(pbm[patchi]))
            {
                const polyPatch& patch = pbm[patchi];
                const vectorField::subField fc = patch.faceCentres();

                // Collect all patch faces
                forAll(fc, patchFacei)
                {
                    label meshFacei = patch.start() + patchFacei;
                    changedFaces[nPatch] = meshFacei;
                    faceDist[nPatch] = wallPoint(fc[patchFacei], 0.0);

                    nPatch++;
                }
            }
        }

        //- Wave calculation engine.
        FaceCellWave<wallPoint> calc
        (
            mesh_,
            changedFaces,
            faceDist,
            allFaceInfo,
            allCellInfo,
            mesh_.globalData().nTotalCells()+1    //maxIter,
        );

        nUnset_ = 0;

        // Copy values
        forAll(allCellInfo, celli)
        {
            scalar dist = allCellInfo[celli].distSqr();

            if (allCellInfo[celli].valid(calc.data()))
            {
                y[celli] = Foam::sqrt(dist);
            }
            else
            {
                y[celli] = dist;
                nUnset_++;
            }
        }
        volScalarField::Boundary& ybf = y.boundaryFieldRef();
        forAll(ybf, patchi)
        {
            fvPatchScalarField& pf = ybf[patchi];
            if (!isA<emptyFvPatchScalarField>(pf))
            {
                const label start = pf.patch().start();
                forAll(pf, patchFacei)
                {
                    label meshFacei = start + patchFacei;
                    scalar dist = allFaceInfo[meshFacei].distSqr();

                    if (allFaceInfo[meshFacei].valid(calc.data()))
                    {
                        // Adding small to avoid problems with /0 in the turbulence
                        // models
                        pf[patchFacei] = Foam::sqrt(dist) + small;
                    }
                    else
                    {
                        pf[patchFacei] = dist;
                        nUnset_++;
                    }
                }
            }
        }
    }
Info<< "** DONE " << endl;

    return nUnset_ > 0;
}


bool Foam::patchDistMethods::segregatedMeshWave::correct
(
    volScalarField& y,
    volVectorField& n
)
{
//     y = dimensionedScalar("yWall", dimLength, great);
// 
//     // Collect pointers to data on patches
//     UPtrList<vectorField> patchData(mesh_.boundaryMesh().size());
// 
//     volVectorField::Boundary& nbf = n.boundaryFieldRef();
// 
//     forAll(nbf, patchi)
//     {
//         patchData.set(patchi, &nbf[patchi]);
//     }
// 
//     // Do mesh wave
//     patchDataWave<wallPointData<vector>> wave
//     (
//         mesh_,
//         patchIDs_,
//         patchData,
//         correctWalls_
//     );
// 
//     // Transfer cell values from wave into y and n
//     y.transfer(wave.distance());
// 
//     n.transfer(wave.cellData());
// 
//     // Transfer values on patches into boundaryField of y and n
//     volScalarField::Boundary& ybf = y.boundaryFieldRef();
// 
//     forAll(ybf, patchi)
//     {
//         scalarField& waveFld = wave.patchDistance()[patchi];
// 
//         if (!isA<emptyFvPatchScalarField>(ybf[patchi]))
//         {
//             ybf[patchi].transfer(waveFld);
// 
//             vectorField& wavePatchData = wave.patchData()[patchi];
// 
//             nbf[patchi].transfer(wavePatchData);
//         }
//     }
// 
//     // Transfer number of unset values
//     nUnset_ = wave.nUnset();
// 
//     return nUnset_ > 0;
    NotImplemented;
    return false;
}


// ************************************************************************* //
