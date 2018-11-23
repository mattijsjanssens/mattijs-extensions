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

#include "procWavePatchDistMethod.H"
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
    defineTypeNameAndDebug(procWave, 0);
    addToRunTimeSelectionTable(patchDistMethod, procWave, dictionary);
}
}

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::patchDistMethods::procWave::furthest
(
    const boundBox& bb,
    const point& pt
) const
{
    point far;
    for (direction dir = 0; dir < pTraits<point>::nComponents; dir++)
    {
        const scalar s = pt[dir];
        const scalar min = bb.min()[dir];
        const scalar max = bb.max()[dir];

        far[dir] = (mag(s-min) > mag(s-max) ? min : max);
    }
    return far;
}


bool Foam::patchDistMethods::procWave::isCloser
(
    const boundBox& localBb,
    const point& localNearest,
    const point& otherNearest
) const
{
    // Get point on bb furthest away
    const point localFar(furthest(localBb, localNearest));
    const point otherFar(furthest(localBb, otherNearest));

    return (magSqr(otherFar-otherNearest) < magSqr(localFar-localNearest));
}


void Foam::patchDistMethods::procWave::processorWave
(
    point& nearestWallPoint
) const
{
    const labelList& neighbourProcs = mesh_.globalData()[Pstream::myProcNo()];
    const boundBox localBb(mesh_.points(), false);

    pointField sendBufs(neighbourProcs.size());
    pointField recvBufs(neighbourProcs.size());

    Pout<< "** STARTING:" << nearestWallPoint << endl;


    while (true)
    {
        // Send my data across
        {
            label startOfRequests = Pstream::nRequests();
            forAll(neighbourProcs, i)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    neighbourProcs[i],
                    reinterpret_cast<char*>(&recvBufs[i]),
                    sizeof(point)
                );
            }
            sendBufs = nearestWallPoint;
            forAll(neighbourProcs, i)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    neighbourProcs[i],
                    reinterpret_cast<const char*>(&sendBufs[i]),
                    sizeof(point)
                );
            }
            Pstream::waitRequests(startOfRequests);
        }

        Pout<< "My data:" << nearestWallPoint << nl
            << "Received:" << recvBufs << endl;

        // Check if any processor is nearer
        bool changed = false;
        forAll(recvBufs, i)
        {
            // Check if wallPoint from neighbour processor is guaranteed nearer
            if (isCloser(localBb, nearestWallPoint, recvBufs[i]))
            {
                Pout<< "Data from processor:" << neighbourProcs[i]
                    << " point:" << recvBufs[i]
                    << " is closer than:" << nearestWallPoint
                    << endl;
                nearestWallPoint = recvBufs[i];
                changed = true;
            }
        }

        if (!returnReduce(changed, orOp<bool>()))
        {
            break;
        }
    }


    Pout<< "** FINAL:" << nearestWallPoint << endl;

    if (debug)
    {
        pointField allNearest(Pstream::nProcs());
        allNearest[Pstream::myProcNo()] = nearestWallPoint;
        Pstream::gatherList(allNearest);

        List<boundBox> allBb(Pstream::nProcs());
        allBb[Pstream::myProcNo()] = localBb;
        Pstream::gatherList(allBb);
        if (Pstream::master())
        {
            OBJstream str("nearest.obj");
            forAll(allBb, proci)
            {
                str.write
                (
                    linePointRef(allBb[proci].midpoint(),
                    allNearest[proci])
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::procWave::procWave
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


Foam::patchDistMethods::procWave::procWave
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

bool Foam::patchDistMethods::procWave::correct(volScalarField& y)
{
    y = dimensionedScalar("yWall", dimLength, great);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    label nPatch = 0;
    forAll(pbm, patchi)
    {
        if (patchIDs_.found(patchi))
        {
            nPatch += pbm[patchi].size();
        }
    }

    // Collect patch faces a (one of) its points

    List<wallPoint> faceDist(nPatch);
    labelList changedFaces(nPatch);

    point nearestWallPoint(vector::uniform(GREAT));

    nPatch = 0;
    forAll(pbm, patchi)
    {
        if (patchIDs_.found(patchi))
        {
            const polyPatch& patch = pbm[patchi];
            const vectorField::subField fc = patch.faceCentres();

            if (fc.size())
            {
                // Pick any face as parallel seed
                nearestWallPoint = fc[0];

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
    }

    //- Wall information for all faces
    List<wallPoint> allFaceInfo(mesh_.nFaces());

    //- Wall information for all cells
    List<wallPoint> allCellInfo(mesh_.nCells());

    //- Determine nearest from any other processors
    processorWave(nearestWallPoint);

    const pointField& faceCentres = mesh_.faceCentres();
    forAll(faceCentres, facei)
    {
        const point& fc = faceCentres[facei];
        allFaceInfo[facei] = wallPoint
        (
            nearestWallPoint,
            magSqr(fc-nearestWallPoint)
        );
    }

    const pointField& cellCentres = mesh_.cellCentres();
    forAll(cellCentres, celli)
    {
        const point& cc = cellCentres[celli];
        allCellInfo[celli] = wallPoint
        (
            nearestWallPoint,
            magSqr(cc-nearestWallPoint)
        );
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

    return nUnset_ > 0;
}


bool Foam::patchDistMethods::procWave::correct
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
