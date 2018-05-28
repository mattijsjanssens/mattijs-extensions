/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "tetMatcher.H"
#include "isoSurface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::isoSurfaceCell::generatePoint
(
    const DynamicList<Type>& snappedPoints,

    const scalar s0,
    const Type& p0,
    const label p0Index,

    const scalar s1,
    const Type& p1,
    const label p1Index
) const
{
    scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        scalar s = (iso_-s0)/d;

        if (s >= 0.5 && s <= 1 && p1Index != -1)
        {
            return snappedPoints[p1Index];
        }
        else if (s >= 0.0 && s <= 0.5 && p0Index != -1)
        {
            return snappedPoints[p0Index];
        }
        else
        {
            return s*p1 + (1.0-s)*p0;
        }
    }
    else
    {
        scalar s = 0.4999;

        return s*p1 + (1.0-s)*p0;
    }
}


//template<class Type>
//void Foam::isoSurfaceCell::generateTriPoints
//(
//    const DynamicList<Type>& snapped,
//
//    const scalar s0,
//    const Type& p0,
//    const label p0Index,
//    const bool isCc0,
//
//    const scalar s1,
//    const Type& p1,
//    const label p1Index,
//    const bool isCc1,
//
//    const scalar s2,
//    const Type& p2,
//    const label p2Index,
//    const bool isCc2,
//
//    const scalar s3,
//    const Type& p3,
//    const label p3Index,
//    const bool isCc3,
//
//    DynamicList<Type>& pts,
//    DynamicList<bool>& usesCc
//) const
//{
//    int triIndex = 0;
//    if (s0 < iso_)
//    {
//        triIndex |= 1;
//    }
//    if (s1 < iso_)
//    {
//        triIndex |= 2;
//    }
//    if (s2 < iso_)
//    {
//        triIndex |= 4;
//    }
//    if (s3 < iso_)
//    {
//        triIndex |= 8;
//    }
//
//    /* Form the vertices of the triangles for each case */
//    switch (triIndex)
//    {
//        case 0x00:
//        case 0x0F:
//        break;
//
//        case 0x01:
//        case 0x0E:
//        {
//            pts.append(generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index));
//            usesCc.append(isCc0||isCc1);
//            pts.append(generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index));
//            usesCc.append(isCc0||isCc2);
//            pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
//            usesCc.append(isCc0||isCc3);
//
//            if (triIndex == 0x0E)
//            {
//                // Flip normals
//                label sz = pts.size();
//                Swap(pts[sz-2], pts[sz-1]);
//                Swap(usesCc[sz-2], usesCc[sz-1]);
//            }
//        }
//        break;
//
//        case 0x02:
//        case 0x0D:
//        {
//            pts.append(generatePoint(snapped,s1,p1,p1Index,s0,p0,p0Index));
//            usesCc.append(isCc1||isCc0);
//            pts.append(generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index));
//            usesCc.append(isCc1||isCc3);
//            pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
//            usesCc.append(isCc1||isCc2);
//
//            if (triIndex == 0x0D)
//            {
//                // Flip normals
//                label sz = pts.size();
//                Swap(pts[sz-2], pts[sz-1]);
//                Swap(usesCc[sz-2], usesCc[sz-1]);
//            }
//        }
//        break;
//
//        case 0x03:
//        case 0x0C:
//        {
//            Type p0p2 = generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index);
//            Type p1p3 = generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index);
//
//            pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
//            usesCc.append(isCc0||isCc3);
//            pts.append(p1p3);
//            usesCc.append(isCc1||isCc3);
//            pts.append(p0p2);
//            usesCc.append(isCc0||isCc2);
//
//            pts.append(p1p3);
//            usesCc.append(isCc1||isCc3);
//            pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
//            usesCc.append(isCc1||isCc2);
//            pts.append(p0p2);
//            usesCc.append(isCc0||isCc2);
//
//            if (triIndex == 0x0C)
//            {
//                // Flip normals
//                label sz = pts.size();
//                Swap(pts[sz-5], pts[sz-4]);
//                Swap(usesCc[sz-5], usesCc[sz-4]);
//                Swap(pts[sz-2], pts[sz-1]);
//                Swap(usesCc[sz-2], usesCc[sz-1]);
//            }
//        }
//        break;
//
//        case 0x04:
//        case 0x0B:
//        {
//            pts.append(generatePoint(snapped,s2,p2,p2Index,s0,p0,p0Index));
//            usesCc.append(isCc2||isCc0);
//            pts.append(generatePoint(snapped,s2,p2,p2Index,s1,p1,p1Index));
//            usesCc.append(isCc2||isCc1);
//            pts.append(generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index));
//            usesCc.append(isCc2||isCc3);
//
//            if (triIndex == 0x0B)
//            {
//                // Flip normals
//                label sz = pts.size();
//                Swap(pts[sz-2], pts[sz-1]);
//                Swap(usesCc[sz-2], usesCc[sz-1]);
//            }
//        }
//        break;
//
//        case 0x05:
//        case 0x0A:
//        {
//            Type p0p1 = generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index);
//            Type p2p3 = generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index);
//
//            pts.append(p0p1);
//            usesCc.append(isCc0||isCc1);
//            pts.append(p2p3);
//            usesCc.append(isCc2||isCc3);
//            pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
//            usesCc.append(isCc0||isCc3);
//
//            pts.append(p0p1);
//            usesCc.append(isCc0||isCc1);
//            pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
//            usesCc.append(isCc1||isCc2);
//            pts.append(p2p3);
//            usesCc.append(isCc2||isCc3);
//
//            if (triIndex == 0x0A)
//            {
//                // Flip normals
//                label sz = pts.size();
//                Swap(pts[sz-5], pts[sz-4]);
//                Swap(usesCc[sz-5], usesCc[sz-4]);
//                Swap(pts[sz-2], pts[sz-1]);
//                Swap(usesCc[sz-2], usesCc[sz-1]);
//            }
//        }
//        break;
//
//        case 0x06:
//        case 0x09:
//        {
//            Type p0p1 = generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index);
//            Type p2p3 = generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index);
//
//            pts.append(p0p1);
//            usesCc.append(isCc0||isCc1);
//            pts.append(generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index));
//            usesCc.append(isCc1||isCc3);
//            pts.append(p2p3);
//            usesCc.append(isCc2||isCc3);
//
//            pts.append(p0p1);
//            usesCc.append(isCc0||isCc1);
//            pts.append(p2p3);
//            usesCc.append(isCc2||isCc3);
//            pts.append(generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index));
//            usesCc.append(isCc0||isCc2);
//
//            if (triIndex == 0x09)
//            {
//                // Flip normals
//                label sz = pts.size();
//                Swap(pts[sz-5], pts[sz-4]);
//                Swap(usesCc[sz-5], usesCc[sz-4]);
//                Swap(pts[sz-2], pts[sz-1]);
//                Swap(usesCc[sz-2], usesCc[sz-1]);
//            }
//        }
//        break;
//
//        case 0x08:
//        case 0x07:
//        {
//            pts.append(generatePoint(snapped,s3,p3,p3Index,s0,p0,p0Index));
//            usesCc.append(isCc3||isCc0);
//            pts.append(generatePoint(snapped,s3,p3,p3Index,s2,p2,p2Index));
//            usesCc.append(isCc3||isCc2);
//            pts.append(generatePoint(snapped,s3,p3,p3Index,s1,p1,p1Index));
//            usesCc.append(isCc3||isCc1);
//            if (triIndex == 0x07)
//            {
//                // Flip normals
//                label sz = pts.size();
//                Swap(pts[sz-2], pts[sz-1]);
//                Swap(usesCc[sz-2], usesCc[sz-1]);
//            }
//        }
//        break;
//    }
//}
//XXXXXXXXXXXX
template<class Type>
void Foam::isoSurfaceCell::generateTriPoints
(
    const DynamicList<Type>& snapped,

    const FixedList<scalar, 4>& s,
    const FixedList<Type, 4>& p,
    const FixedList<label, 4>& pIndex,

    const label ccIndex,                 // index (tet numbering) of cc
    const FixedList<bool, 6>& edgeIsDiag,// per tet edge whether is face diag

    DynamicList<Type>& pts,
    DynamicList<bool>& usesCc,
    DynamicList<bool>& onDiag
) const
{
    int triIndex = 0;
    if (s[0] < iso_)
    {
        triIndex |= 1;
    }
    if (s[1] < iso_)
    {
        triIndex |= 2;
    }
    if (s[2] < iso_)
    {
        triIndex |= 4;
    }
    if (s[3] < iso_)
    {
        triIndex |= 8;
    }

    /* Form the vertices of the triangles for each case */
    switch (triIndex)
    {
        case 0x00:
        case 0x0F:
        break;

        case 0x01:
        case 0x0E:
        {
            pts.append
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[1],p[1],pIndex[1])
            );
            usesCc.append(ccIndex==0||ccIndex==1);
            onDiag.append(edgeIsDiag[0]);
            pts.append
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[2],p[2],pIndex[2])
            );
            usesCc.append(ccIndex==0||ccIndex==2);
            onDiag.append(edgeIsDiag[1]);
            pts.append
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[3],p[3],pIndex[3])
            );
            usesCc.append(ccIndex==0||ccIndex==3);
            onDiag.append(edgeIsDiag[2]);

            if (triIndex == 0x0E)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-2], pts[sz-1]);
                Swap(usesCc[sz-2], usesCc[sz-1]);
                Swap(onDiag[sz-2], onDiag[sz-1]);
            }
        }
        break;

        case 0x02:
        case 0x0D:
        {
            pts.append
            (
                generatePoint(snapped,s[1],p[1],pIndex[1],s[0],p[0],pIndex[0])
            );
            usesCc.append(ccIndex==1||ccIndex==0);
            onDiag.append(edgeIsDiag[0]);
            pts.append
            (
                generatePoint(snapped,s[1],p[1],pIndex[1],s[3],p[3],pIndex[3])
            );
            usesCc.append(ccIndex==1||ccIndex==3);
            onDiag.append(edgeIsDiag[3]);
            pts.append
            (
                generatePoint(snapped,s[1],p[1],pIndex[1],s[2],p[2],pIndex[2])
            );
            usesCc.append(ccIndex==1||ccIndex==2);
            onDiag.append(edgeIsDiag[4]);

            if (triIndex == 0x0D)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-2], pts[sz-1]);
                Swap(usesCc[sz-2], usesCc[sz-1]);
                Swap(onDiag[sz-2], onDiag[sz-1]);
            }
        }
        break;

        case 0x03:
        case 0x0C:
        {
            Type p0p2
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[2],p[2],pIndex[2])
            );
            Type p1p3
            (
                generatePoint(snapped,s[1],p[1],pIndex[1],s[3],p[3],pIndex[3])
            );

            pts.append
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[3],p[3],pIndex[3])
            );
            usesCc.append(ccIndex==0||ccIndex==3);
            onDiag.append(edgeIsDiag[2]);
            pts.append(p1p3);
            usesCc.append(ccIndex==1||ccIndex==3);
            onDiag.append(edgeIsDiag[3]);
            pts.append(p0p2);
            usesCc.append(ccIndex==0||ccIndex==2);
            onDiag.append(edgeIsDiag[1]);

            pts.append(p1p3);
            usesCc.append(ccIndex==1||ccIndex==3);
            onDiag.append(edgeIsDiag[3]);
            pts.append
            (
                generatePoint(snapped,s[1],p[1],pIndex[1],s[2],p[2],pIndex[2])
            );
            usesCc.append(ccIndex==1||ccIndex==2);
            onDiag.append(edgeIsDiag[4]);
            pts.append(p0p2);
            usesCc.append(ccIndex==0||ccIndex==2);
            onDiag.append(edgeIsDiag[1]);

            if (triIndex == 0x0C)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-5], pts[sz-4]);
                Swap(usesCc[sz-5], usesCc[sz-4]);
                Swap(onDiag[sz-5], onDiag[sz-4]);
                Swap(pts[sz-2], pts[sz-1]);
                Swap(usesCc[sz-2], usesCc[sz-1]);
                Swap(onDiag[sz-2], onDiag[sz-1]);
            }
        }
        break;

        case 0x04:
        case 0x0B:
        {
            pts.append
            (
                generatePoint(snapped,s[2],p[2],pIndex[2],s[0],p[0],pIndex[0])
            );
            usesCc.append(ccIndex==2||ccIndex==0);
            onDiag.append(edgeIsDiag[1]);
            pts.append
            (
                generatePoint(snapped,s[2],p[2],pIndex[2],s[1],p[1],pIndex[1])
            );
            usesCc.append(ccIndex==2||ccIndex==1);
            onDiag.append(edgeIsDiag[4]);
            pts.append
            (
                generatePoint(snapped,s[2],p[2],pIndex[2],s[3],p[3],pIndex[3])
            );
            usesCc.append(ccIndex==2||ccIndex==3);
            onDiag.append(edgeIsDiag[5]);

            if (triIndex == 0x0B)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-2], pts[sz-1]);
                Swap(usesCc[sz-2], usesCc[sz-1]);
                Swap(onDiag[sz-2], onDiag[sz-1]);
            }
        }
        break;

        case 0x05:
        case 0x0A:
        {
            Type p0p1
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[1],p[1],pIndex[1])
            );
            Type p2p3
            (
                generatePoint(snapped,s[2],p[2],pIndex[2],s[3],p[3],pIndex[3])
            );

            pts.append(p0p1);
            usesCc.append(ccIndex==0||ccIndex==1);
            onDiag.append(edgeIsDiag[0]);
            pts.append(p2p3);
            usesCc.append(ccIndex==2||ccIndex==3);
            onDiag.append(edgeIsDiag[5]);
            pts.append
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[3],p[3],pIndex[3]))
            );
            usesCc.append(ccIndex==0||ccIndex==3);
            onDiag.append(edgeIsDiag[2]);

            pts.append(p0p1);
            usesCc.append(ccIndex==0||ccIndex==1);
            onDiag.append(edgeIsDiag[0]);
            pts.append
            (
                generatePoint(snapped,s[1],p[1],pIndex[1],s[2],p[2],pIndex[2])
            );
            usesCc.append(ccIndex==1||ccIndex==2);
            onDiag.append(edgeIsDiag[4]);
            pts.append(p2p3);
            usesCc.append(ccIndex==2||ccIndex==3);
            onDiag.append(edgeIsDiag[5]);

            if (triIndex == 0x0A)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-5], pts[sz-4]);
                Swap(usesCc[sz-5], usesCc[sz-4]);
                Swap(onDiag[sz-5], onDiag[sz-4]);
                Swap(pts[sz-2], pts[sz-1]);
                Swap(usesCc[sz-2], usesCc[sz-1]);
                Swap(onDiag[sz-2], onDiag[sz-1]);
            }
        }
        break;

        case 0x06:
        case 0x09:
        {
            Type p0p1
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[1],p[1],pIndex[1])
            );
            Type p2p3
            (
                generatePoint(snapped,s[2],p[2],pIndex[2],s[3],p[3],pIndex[3])
            );

            pts.append(p0p1);
            usesCc.append(ccIndex==0||ccIndex==1);
            onDiag.append(edgeIsDiag[0]);
            pts.append
            (
                generatePoint(snapped,s[1],p[1],pIndex[1],s[3],p[3],pIndex[3])
            );
            usesCc.append(ccIndex==1||ccIndex==3);
            onDiag.append(edgeIsDiag[3]);
            pts.append(p2p3);
            usesCc.append(ccIndex==2||ccIndex==3);
            onDiag.append(edgeIsDiag[5]);

            pts.append(p0p1);
            usesCc.append(ccIndex==0||ccIndex==1);
            onDiag.append(edgeIsDiag[0]);
            pts.append(p2p3);
            usesCc.append(ccIndex==2||ccIndex==3);
            onDiag.append(edgeIsDiag[5]);
            pts.append
            (
                generatePoint(snapped,s[0],p[0],pIndex[0],s[2],p[2],pIndex[2])
            );
            usesCc.append(ccIndex==0||ccIndex==2);
            onDiag.append(edgeIsDiag[1]);

            if (triIndex == 0x09)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-5], pts[sz-4]);
                Swap(usesCc[sz-5], usesCc[sz-4]);
                Swap(onDiag[sz-5], onDiag[sz-4]);
                Swap(pts[sz-2], pts[sz-1]);
                Swap(usesCc[sz-2], usesCc[sz-1]);
                Swap(onDiag[sz-2], onDiag[sz-1]);
            }
        }
        break;

        case 0x08:
        case 0x07:
        {
            pts.append
            (
                generatePoint(snapped,s[3],p[3],pIndex[3],s[0],p[0],pIndex[0])
            );
            usesCc.append(ccIndex==3||ccIndex==0);
            onDiag.append(edgeIsDiag[2]);
            pts.append
            (
                generatePoint(snapped,s[3],p[3],pIndex[3],s[2],p[2],pIndex[2])
            );
            usesCc.append(ccIndex==3||ccIndex==2);
            onDiag.append(edgeIsDiag[5]);
            pts.append
            (
                generatePoint(snapped,s[3],p[3],pIndex[3],s[1],p[1],pIndex[1])
            );
            usesCc.append(ccIndex==3||ccIndex==1);
            onDiag.append(edgeIsDiag[3]);
            if (triIndex == 0x07)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-2], pts[sz-1]);
                Swap(usesCc[sz-2], usesCc[sz-1]);
                Swap(onDiag[sz-2], onDiag[sz-1]);
            }
        }
        break;
    }
}
//XXXXXXX


template<class Type>
void Foam::isoSurfaceCell::generateTriPoints
(
    const scalarField& cVals,
    const scalarField& pVals,

    const Field<Type>& cCoords,
    const Field<Type>& pCoords,

    const DynamicList<Type>& snappedPoints,
    const labelList& snappedCc,
    const labelList& snappedPoint,

    DynamicList<Type>& triPoints,
    DynamicList<label>& triMeshCells,
    DynamicList<bool>& usesCellCentre,
    DynamicList<bool>& usesDiag
) const
{
    tetMatcher tet;
    label countNotFoundTets = 0;

    FixedList<scalar, 4> s;
    FixedList<Type, 4> p;
    FixedList<label, 4> pIndex;
    FixedList<bool, 6> edgeIsDiag;

    forAll(mesh_.cells(), celli)
    {
        if (cellCutType_[celli] != NOTCUT)
        {
            label oldNPoints = triPoints.size();

            const cell& cFaces = mesh_.cells()[celli];

            if (tet.isA(mesh_, celli))
            {
                // For tets don't do cell-centre decomposition, just use the
                // tet points and values

                const face& f0 = mesh_.faces()[cFaces[0]];

                // Get the other point
                const face& f1 = mesh_.faces()[cFaces[1]];
                label oppositeI = -1;
                forAll(f1, fp)
                {
                    oppositeI = f1[fp];

                    if (findIndex(f0, oppositeI) == -1)
                    {
                        break;
                    }
                }

                // Start off from positive volume tet to make sure we
                // generate outwards pointing tets
                if (mesh_.faceOwner()[cFaces[0]] == celli)
                {
                    s[0] = pVals[f0[1]];
                    p[0] = pCoords[f0[1]];
                    pIndex[0] = snappedPoint[f0[1]];

                    s[1] = pVals[f0[0]];
                    p[1] = pCoords[f0[0]];
                    pIndex[1] = snappedPoint[f0[0]];

                    s[2] = pVals[f0[2]];
                    p[2] = pCoords[f0[2]];
                    pIndex[2] = snappedPoint[f0[2]];

                    s[3] = pVals[oppositeI];
                    p[3] = pCoords[oppositeI];
                    pIndex[3] = snappedPoint[oppositeI];

                    edgeIsDiag = false;

                    generateTriPoints
                    (
                        snappedPoints,

                        s,
                        p,
                        pIndex,

                        -1,             // index in tet numbering of cell centre
                        edgeIsDiag,     // per tet edge whether is face diag

                        triPoints,
                        usesCellCentre,
                        usesDiag
                    );
                }
                else
                {
                    s[0] = pVals[f0[0]];
                    p[0] = pCoords[f0[0]];
                    pIndex[0] = snappedPoint[f0[0]];

                    s[1] = pVals[f0[1]];
                    p[1] = pCoords[f0[1]];
                    pIndex[1] = snappedPoint[f0[1]];

                    s[2] = pVals[f0[2]];
                    p[2] = pCoords[f0[2]];
                    pIndex[2] = snappedPoint[f0[2]];

                    s[3] = pVals[oppositeI];
                    p[3] = pCoords[oppositeI];
                    pIndex[3] = snappedPoint[oppositeI];

                    edgeIsDiag = false;

                    generateTriPoints
                    (
                        snappedPoints,

                        s,
                        p,
                        pIndex,

                        -1,             // index in tet numbering of cell centre
                        edgeIsDiag,     // per tet edge whether is face diag

                        triPoints,
                        usesCellCentre,
                        usesDiag
                    );
                }
            }
            else
            {
                forAll(cFaces, cFacei)
                {
                    label facei = cFaces[cFacei];
                    const face& f = mesh_.faces()[facei];

                    label fp0 = mesh_.tetBasePtIs()[facei];

                    // Skip undefined tets
                    if (fp0 < 0)
                    {
                        fp0 = 0;
                        countNotFoundTets++;
                    }

                    label fp = f.fcIndex(fp0);
                    for (label i = 2; i < f.size(); i++)
                    {
                        label nextFp = f.fcIndex(fp);
                        triFace tri(f[fp0], f[fp], f[nextFp]);

                        // Start off from positive volume tet to make sure we
                        // generate outwards pointing tets
                        if (mesh_.faceOwner()[facei] == celli)
                        {
                            s[0] = pVals[tri[1]];
                            p[0] = pCoords[tri[1]];
                            pIndex[0] = snappedPoint[tri[1]];

                            s[1] = pVals[tri[0]];
                            p[1] = pCoords[tri[0]];
                            pIndex[1] = snappedPoint[tri[0]];

                            s[2] = pVals[tri[2]];
                            p[2] = pCoords[tri[2]];
                            pIndex[2] = snappedPoint[tri[2]];

                            s[3] = cVals[celli];
                            p[3] = cCoords[celli];
                            pIndex[3] = snappedCc[celli];

                            edgeIsDiag = false;
                            if (i != 2) edgeIsDiag[0] = true;
                            if (i != f.size()-1) edgeIsDiag[4] = true;

                            generateTriPoints
                            (
                                snappedPoints,

                                s,
                                p,
                                pIndex,

                                3,          // index of cell centre
                                edgeIsDiag, // per tet edge whether is face diag

                                triPoints,
                                usesCellCentre,
                                usesDiag
                            );
                        }
                        else
                        {
                            s[0] = pVals[tri[0]];
                            p[0] = pCoords[tri[0]];
                            pIndex[0] = snappedPoint[tri[0]];

                            s[1] = pVals[tri[1]];
                            p[1] = pCoords[tri[1]];
                            pIndex[1] = snappedPoint[tri[1]];

                            s[2] = pVals[tri[2]];
                            p[2] = pCoords[tri[2]];
                            pIndex[2] = snappedPoint[tri[2]];

                            s[3] = cVals[celli];
                            p[3] = cCoords[celli];
                            pIndex[3] = snappedCc[celli];

                            edgeIsDiag = false;
                            if (i != 2) edgeIsDiag[0] = true;
                            if (i != f.size()-1) edgeIsDiag[1] = true;

                            generateTriPoints
                            (
                                snappedPoints,

                                s,
                                p,
                                pIndex,

                                3,          // cell centre index
                                edgeIsDiag, // per tet edge :is face diag

                                triPoints,
                                usesCellCentre,
                                usesDiag
                            );
                        }

                        fp = nextFp;
                    }
                }
            }

            // Every three triPoints is a cell
            label nCells = (triPoints.size()-oldNPoints)/3;
            for (label i = 0; i < nCells; i++)
            {
                triMeshCells.append(celli);
            }
        }
    }

    if (countNotFoundTets > 0)
    {
        WarningInFunction
            << "Could not find " << countNotFoundTets
            << " tet base points, which may lead to inverted triangles."
            << endl;
    }

    triPoints.shrink();
    triMeshCells.shrink();
    usesCellCentre.shrink();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurfaceCell::interpolate
(
    const Field<Type>& cCoords,
    const Field<Type>& pCoords
) const
{
    DynamicList<Type> triPoints(3*nCutCells_);
    DynamicList<label> triMeshCells(nCutCells_);
    DynamicList<bool> usesCellCentre(3*nCutCells_);
    DynamicList<bool> usesDiag(3*nCutCells_);

    // Dummy snap data
    DynamicList<Type> snappedPoints;
    labelList snappedCc(mesh_.nCells(), -1);
    labelList snappedPoint(mesh_.nPoints(), -1);

    generateTriPoints
    (
        cVals_,
        pVals_,

        cCoords,
        pCoords,

        snappedPoints,
        snappedCc,
        snappedPoint,

        triPoints,
        triMeshCells,
        usesCellCentre,
        usesDiag
    );

    return isoSurface::interpolate
    (
        points().size(),
        triPointMergeMap_,
        triPoints
    );
}


// ************************************************************************* //
