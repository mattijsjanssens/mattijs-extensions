/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "isoSurfaceCellTopo.H"
// #include "dictionary.H"
// #include "polyMesh.H"
// #include "mergePoints.H"
// #include "tetMatcher.H"
// #include "syncTools.H"
// #include "addToRunTimeSelectionTable.H"
// #include "DynamicField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoSurfaceCellTopo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::isoSurfaceCellTopo::isoFraction
(
    const scalar s0,
    const scalar s1
) const
{
    scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        return (iso_-s0)/d;
    }
    else
    {
        return -1.0;
    }
}


bool Foam::isoSurfaceCellTopo::isTriCut
(
    const triFace& tri,
    const scalarField& pointValues
) const
{
    bool aLower = (pointValues[tri[0]] < iso_);
    bool bLower = (pointValues[tri[1]] < iso_);
    bool cLower = (pointValues[tri[2]] < iso_);

    return !(aLower == bLower && aLower == cLower);
}


Foam::isoSurfaceCellTopo::cellCutType Foam::isoSurfaceCellTopo::calcCutType
(
    const bool isTet,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];

    if (isTet)
    {
        forAll(cFaces, cFacei)
        {
            const face& f = mesh_.faces()[cFaces[cFacei]];

            for (label fp = 1; fp < f.size() - 1; fp++)
            {
                triFace tri(f[0], f[fp], f[f.fcIndex(fp)]);

                if (isTriCut(tri, pVals_))
                {
                    return CUT;
                }
            }
        }
        return NOTCUT;
    }
    else
    {
        bool cellLower = (cVals_[celli] < iso_);

        // First check if there is any cut in cell
        bool edgeCut = false;

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const face& f = mesh_.faces()[facei];

            // Check pyramids cut
            forAll(f, fp)
            {
                if ((pVals_[f[fp]] < iso_) != cellLower)
                {
                    edgeCut = true;
                    break;
                }
            }

            if (edgeCut)
            {
                break;
            }

            const label fp0 = mesh_.tetBasePtIs()[facei];
            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                label nextFp = f.fcIndex(fp);

                if (isTriCut(triFace(f[fp0], f[fp], f[nextFp]), pVals_))
                {
                    edgeCut = true;
                    break;
                }

                fp = nextFp;
            }

            if (edgeCut)
            {
                break;
            }
        }

        if (edgeCut)
        {
            // Count actual cuts (expensive since addressing needed)
            // Note: not needed if you don't want to preserve maxima/minima
            // centred around cellcentre. In that case just always return CUT

            const labelList& cPoints = mesh_.cellPoints(celli);

            label nPyrCuts = 0;

            forAll(cPoints, i)
            {
                if ((pVals_[cPoints[i]] < iso_) != cellLower)
                {
                    nPyrCuts++;
                }
            }

            if (nPyrCuts == cPoints.size())
            {
                return SPHERE;
            }
            else
            {
                return CUT;
            }
        }
        else
        {
            return NOTCUT;
        }
    }
}


void Foam::isoSurfaceCellTopo::calcCutTypes(const tetMatcher& tet)
{
    cellCutType_.setSize(mesh_.nCells());
    nCutCells_ = 0;
    forAll(mesh_.cells(), celli)
    {
        cellCutType_[celli] = calcCutType
        (
            tet.isA(mesh_, celli),
            celli
        );

        if (cellCutType_[celli] == CUT)
        {
            nCutCells_++;
        }
    }

    if (debug)
    {
        Pout<< "isoSurfaceCellTopo : detected " << nCutCells_
            << " candidate cut cells." << endl;
    }
}


Foam::label Foam::isoSurfaceCellTopo::generatePoint
(
    const label facei,
    const bool edgeIsDiag,
    const edge& vertices,

    DynamicList<edge>& pointToVerts,
    DynamicList<label>& pointToFace,
    DynamicList<bool>& pointFromDiag,
    EdgeMap<label>& vertsToPoint
) const
{
    EdgeMap<label>::const_iterator edgeFnd = vertsToPoint.find(vertices);
    if (edgeFnd != vertsToPoint.end())
    {
        return edgeFnd();
    }
    else
    {
        // Generate new point
        label pointi = pointToVerts.size();

        pointToVerts.append(vertices);
        pointToFace.append(facei);
        pointFromDiag.append(edgeIsDiag);
        vertsToPoint.insert(vertices, pointi);
        return pointi;
    }
}


void Foam::isoSurfaceCellTopo::generateTriPoints
(
    const label facei,
    const FixedList<scalar, 4>& s,
    const FixedList<point, 4>& p,
    const FixedList<label, 4>& pIndex,
    const FixedList<bool, 6>& edgeIsDiag,// per tet edge whether is face diag

    DynamicList<edge>& pointToVerts,
    DynamicList<label>& pointToFace,
    DynamicList<bool>& pointFromDiag,

    EdgeMap<label>& vertsToPoint,
    DynamicList<label>& verts       // every three verts is new triangle
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
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    edge(pIndex[0], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    edge(pIndex[0], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    edge(pIndex[0], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (triIndex == 0x0E)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x02:
        case 0x0D:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    edge(pIndex[1], pIndex[0]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    edge(pIndex[1], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    edge(pIndex[1], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (triIndex == 0x0D)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x03:
        case 0x0C:
        {
            label p0p2
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    edge(pIndex[0], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            label p1p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    edge(pIndex[1], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    edge(pIndex[0], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p1p3);
            verts.append(p0p2);
            verts.append(p1p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    edge(pIndex[1], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p0p2);

            if (triIndex == 0x0C)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-5], verts[sz-4]);
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x04:
        case 0x0B:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    edge(pIndex[2], pIndex[0]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    edge(pIndex[2], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    edge(pIndex[2], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (triIndex == 0x0B)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x05:
        case 0x0A:
        {
            label p0p1
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    edge(pIndex[0], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            label p2p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    edge(pIndex[2], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append(p0p1);
            verts.append(p2p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    edge(pIndex[0], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p0p1);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[4],
                    edge(pIndex[1], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p2p3);

            if (triIndex == 0x0A)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-5], verts[sz-4]);
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x06:
        case 0x09:
        {
            label p0p1
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[0],
                    edge(pIndex[0], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            label p2p3
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    edge(pIndex[2], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            verts.append(p0p1);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    edge(pIndex[1], pIndex[3]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append(p2p3);
            verts.append(p0p1);
            verts.append(p2p3);
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[1],
                    edge(pIndex[0], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );

            if (triIndex == 0x09)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-5], verts[sz-4]);
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;

        case 0x08:
        case 0x07:
        {
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[2],
                    edge(pIndex[3], pIndex[0]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[5],
                    edge(pIndex[3], pIndex[2]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            verts.append
            (
                generatePoint
                (
                    facei,
                    edgeIsDiag[3],
                    edge(pIndex[3], pIndex[1]),
                    pointToVerts, pointToFace, pointFromDiag, vertsToPoint
                )
            );
            if (triIndex == 0x07)
            {
                // Flip normals
                label sz = verts.size();
                Swap(verts[sz-2], verts[sz-1]);
            }
        }
        break;
    }
}


void Foam::isoSurfaceCellTopo::generateTriPoints
(
    const polyMesh& mesh,
    const label celli,
    const bool isTet,

    DynamicList<edge>& pointToVerts,
    DynamicList<label>& pointToFace,
    DynamicList<bool>& pointFromDiag,

    EdgeMap<label>& vertsToPoint,
    DynamicList<label>& verts,
    DynamicList<label>& faceLabels
) const
{
    const cell& cFaces = mesh.cells()[celli];

    if (isTet)
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


        label p0 = f0[0];
        label p1 = f0[1];
        label p2 = f0[2];
        FixedList<bool, 6> edgeIsDiag(false);

        if (mesh.faceOwner()[facei] == celli)
        {
            Swap(p1, p2);
        }

        tetPointRef tet
        (
            mesh.points()[p0],
            mesh.points()[p1],
            mesh.points()[p2],
            mesh.points()[oppositeI]
        );

        label startTrii = verts.size();
        generateTriPoints
        (
            facei,
            FixedList<scalar, 4>
            ({
                pVals_[p0],
                pVals_[p1],
                pVals_[p2],
                pVals_[oppositeI]
            }),
            FixedList<point, 4>
            ({
                mesh.points()[p0],
                mesh.points()[p1],
                mesh.points()[p2],
                mesh.points()[oppositeI]
            }),
            FixedList<label, 4>
            ({
                p0,
                p1,
                p2,
                oppositeI
            }),
            edgeIsDiag,

            pointToVerts,
            pointToFace,
            pointFromDiag,
            vertsToPoint,
            verts       // every three verts is new triangle
        );

        label nTris = (verts.size()-startTrii)/3;
        for (label i = 0; i < nTris; i++)
        {
            faceLabels.append(facei);
        }
    }
    else
    {
        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const face& f = mesh.faces()[facei];

            label fp0 = mesh.tetBasePtIs()[facei];

            label startTrii = verts.size();

            // Skip undefined tets
            if (fp0 < 0)
            {
                fp0 = 0;
            }

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                label nextFp = f.fcIndex(fp);

                FixedList<bool, 6> edgeIsDiag(false);

                label p0 = f[fp0];
                label p1 = f[fp];
                label p2 = f[nextFp];
                if (mesh.faceOwner()[facei] == celli)
                {
                    Swap(p1, p2);
                    if (i != 2) edgeIsDiag[1] = true;
                    if (i != f.size()-1) edgeIsDiag[0] = true;
                }
                else
                {
                    if (i != 2) edgeIsDiag[0] = true;
                    if (i != f.size()-1) edgeIsDiag[1] = true;
                }

                tetPointRef tet
                (
                    mesh.points()[p0],
                    mesh.points()[p1],
                    mesh.points()[p2],
                    mesh.cellCentres()[celli]
                );

                generateTriPoints
                (
                    facei,
                    FixedList<scalar, 4>
                    ({
                        pVals_[p0],
                        pVals_[p1],
                        pVals_[p2],
                        cVals_[celli]
                    }),
                    FixedList<point, 4>
                    ({
                        mesh.points()[p0],
                        mesh.points()[p1],
                        mesh.points()[p2],
                        mesh.cellCentres()[celli]
                    }),
                    FixedList<label, 4>
                    ({
                        p0,
                        p1,
                        p2,
                        mesh.nPoints()+celli
                    }),
                    edgeIsDiag,

                    pointToVerts,
                    pointToFace,
                    pointFromDiag,
                    vertsToPoint,
                    verts       // every three verts is new triangle
                );

                fp = nextFp;
            }

            label nTris = (verts.size()-startTrii)/3;
            for (label i = 0; i < nTris; i++)
            {
                faceLabels.append(facei);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceCellTopo::isoSurfaceCellTopo
(
    const polyMesh& mesh,
    const scalarField& cVals,
    const scalarField& pVals,
    const scalar iso,
    const bool regularise
)
:
    mesh_(mesh),
    cVals_(cVals),
    pVals_(pVals),
    iso_(iso)
{
    if (debug)
    {
        Pout<< "isoSurfaceCellTopo : iso:" << iso_
            << " regularise:" << regularise
            << endl;
    }

    const tetMatcher tet;

    // Determine if any cut through cell
    calcCutTypes(tet, cVals_, pVals_);

    // Per cell: 5 pyramids cut, each generating 2 triangles
    DynamicList<edge> pointToVerts(10*nCutCells_);
    DynamicList<label> pointToFace(10*nCutCells_);
    DynamicList<bool> pointFromDiag(10*nCutCells_);

    // Per cell: number of intersected edges:
    //          - four faces cut so 4 mesh edges + 4 face-diagonal edges
    //          - 4 of the pyramid edges
    EdgeMap<label> vertsToPoint(12*nCutCells_);
    DynamicList<label> verts(12*nCutCells_);
    DynamicList<label> faceLabels(5*nCutCells_);
    DynamicList<label> cellLabels(5*nCutCells_);

    for (label celli = 0; celli < mesh_.nCells(); celli++)
    {
        if (cellCutType_[celli] != NOTCUT)
        {
            label nOldTris = faceLabels.size();

            generateTriPoints
            (
                mesh,
                celli,
                tet.isA(mesh_, celli),

                pointToVerts,
                pointToFace,
                pointFromDiag,

                vertsToPoint,
                verts,
                faceLabels
            );

            for (label i = nOldTris; i < faceLabels.size(); i++)
            {
                cellLabels.append(celli);
            }
        }
    }

    pointToVerts_.transfer(pointToVerts);
    meshCells_.transfer(cellLabels);

    tmp<pointField> allPoints
    (
        interpolate
        (
            mesh_,
            mesh_.points()
        )
    );

    // Assign to MeshedSurface
    faceList allTris(faceLabels.size());
    label verti = 0;
    forAll(allTris, i)
    {
        allTris[i].setSize(3);
        allTris[i][0] = verts[verti++];
        allTris[i][1] = verts[verti++];
        allTris[i][2] = verts[verti++];
    }



//XXXXXXX



        if (debug)
        {
            Pout<< "isoSurfaceCellTopo : generated " << triMeshCells.size()
                << " unmerged triangles." << endl;
        }

        // Merge points and compact out non-valid triangles
        labelList triMap;
        triSurface::operator=
        (
            stitchTriPoints
            (
                false,              // check for duplicate tris
                triPoints,
                triPointMergeMap_,  // unmerged to merged point
                triMap              // merged to unmerged triangle
            )
        );

        if (debug)
        {
            Pout<< "isoSurfaceCellTopo : generated " << triMap.size()
                << " merged triangles." << endl;
        }

        meshCells_.setSize(triMap.size());
        forAll(triMap, i)
        {
            meshCells_[i] = triMeshCells[triMap[i]];
        }
        // Very close non-diag and diag points might get merged. If so
        // make sure non-diag wins
        isOnDiag_.setSize(this->points().size());
        isOnDiag_ = true;
        forAll(triPointMergeMap_, i)
        {
            if (!usesDiag[i])
            {
                isOnDiag_[triPointMergeMap_[i]] = false;
            }
        }
        const boolList oldUsesCellCentre(usesCellCentre.xfer());
        usesCellCentre.setSize(points().size());
        usesCellCentre = true;
        forAll(triPointMergeMap_, i)
        {
            if (!oldUsesCellCentre[i])
            {
                usesCellCentre[triPointMergeMap_[i]] = false;
            }
        }

        if (size() && regularise)
        {
            const label nOldPoints = points().size();

            // Triangulate outside
            DynamicList<label> pointCompactMap; // back to original point
            DynamicList<label> compactCellIDs;  // per returned tri the cellID
            triSurface::operator=
            (
                removeInsidePoints
                (
                    true,   //removeDiagPoints
                    *this,
                    isOnDiag_,
                    meshCells_,
                    pointCompactMap,
                    compactCellIDs
                )
            );

            if (debug)
            {
                Pout<< "isoSurfaceCellTopo :"
                    << " after removing cell centre triangles : " << size()
                    << endl;
            }

            meshCells_.transfer(compactCellIDs);
            labelList reversePointMap(invert(nOldPoints, pointCompactMap));
            inplaceRenumber(reversePointMap, triPointMergeMap_);
            isOnDiag_ = UIndirectList<bool>(isOnDiag_, pointCompactMap)();
        }
    }


    if (debug)
    {
        Pout<< "isoSurfaceCellTopo : checking " << size()
            << " triangles for validity." << endl;

        forAll(*this, triI)
        {
            // Copied from surfaceCheck
            validTri(*this, triI);
        }
    }
}


// ************************************************************************* //
