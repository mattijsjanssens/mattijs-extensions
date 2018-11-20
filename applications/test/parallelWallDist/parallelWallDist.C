/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 208 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "fvMesh.H"
#include "wallPointData.H"
#include "FaceCellWave.H"
#include "fvcGrad.H"
#include "localMin.H"
#include "localMax.H"
#include "OBJstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//vector maxDistance
//(
//    const boundBox& myBb,
//    const boundBox& otherBb
//)
//{
//    // Check distance in x
//    vector maxDist(cmptMax(myBb.min()-otherBb.min()));
//    maxDist = cmptMax(maxDist, myBb.max()-otherBb.min());
//    maxDist = cmptMax(maxDist, myBb.max()-otherBb.max());
//    maxDist = cmptMax(maxDist, myBb.min()-otherBb.max());
//
//    return maxDist;
//}


point furthest(const boundBox& bb, const point& pt)
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


bool isCloser
(
    const boundBox& localBb,
    const point& localNearest,
    const point& otherNearest
)
{
    // Get point on bb furthest away
    const point localFar(furthest(localBb, localNearest));
    const point otherFar(furthest(localBb, otherNearest));

    return (magSqr(otherFar-otherNearest) < magSqr(localFar-localNearest));
}


int main(int argc, char *argv[])
{
    argList::validArgs.append("patches");

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"

    const boundBox localBb(mesh.points(), false);


    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Find set of patches from the list of regular expressions provided
    const wordReList patches((IStringStream(args[1])()));
    const labelHashSet patchSet(pbm.patchSet(patches));

    if (!patchSet.size())
    {
        FatalErrorInFunction
            << "Cannot find any patches in set " << patches << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    // Check my local size of patches
    bool havePatch = false;
    point nearestWallPoint(vector::uniform(GREAT));
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        const polyPatch& pp = pbm[iter.key()];
        if (pp.size())
        {
            nearestWallPoint = pp.points()[pp[0][0]];
            havePatch = true;
            break;
        }
    }

DebugVar(havePatch);

    const labelList& neighbourProcs = mesh.globalData()[Pstream::myProcNo()];
    pointField sendBufs(neighbourProcs.size());
    pointField recvBufs(neighbourProcs.size());

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
            str.write(linePointRef(allBb[proci].midpoint(), allNearest[proci]));
        }
    }

//
//    volScalarField fld
//    (
//        IOobject
//        (
//            "interfaceHeight",
//            mesh.time().timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE,
//            false
//        ),
//        mesh,
//        dimensionedScalar("large", dimless, mesh.globalData().nTotalCells()+1)
//    );
//    forAll(cellData, celli)
//    {
//        if (cellData[celli].valid(deltaCalc.data()))
//        {
//            fld[celli] = cellData[celli].data();
//        }
//    }
//    forAll(fld.boundaryFieldRef(), patchi)
//    {
//        fvPatchScalarField& fvp = fld.boundaryFieldRef()[patchi];
//        scalarField patchVals(fvp.size(), 0.0);
//
//        forAll(patchVals, i)
//        {
//            label facei = fvp.patch().start()+i;
//            if (faceData[facei].valid(deltaCalc.data()))
//            {
//                patchVals[i] = faceData[facei].data();
//            }
//        }
//
//        fvp == patchVals;
//    }
//
//    Info<< "Writing " << fld.name() << " with interfaceHeight" << endl;
//    fld.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
