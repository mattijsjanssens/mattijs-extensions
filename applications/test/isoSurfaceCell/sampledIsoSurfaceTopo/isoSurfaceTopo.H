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

Class
    Foam::isoSurfaceTopo

Description
    Marching tet iso surface algorithm with optional filtering to keep only
    points originating from mesh edges.

SourceFiles
    isoSurfaceTopo.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceTopo_H
#define isoSurfaceTopo_H

#include "labelPair.H"
#include "pointIndexHit.H"
#include "PackedBoolList.H"
#include "MeshedSurface.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class tetMatcher;

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceTopo Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceTopo
:
    public MeshedSurface<face>
{
public:

        enum filterType
        {
            NONE,       // no filtering
            DIAGCELL,   // remove points from face-diagonal and pyramid
                        // (vertex to cell-centre) edges
            CELL        // only remove points from pyramid edges
        };


private:

    // Private data

        enum cellCutType
        {
            NOTCUT,     // not cut
            SPHERE,     // all edges to cell centre cut
            CUT         // normal cut
        };


        //- Reference to mesh
        const polyMesh& mesh_;

        const scalarField& cVals_;

        const scalarField& pVals_;

        //- Iso value
        const scalar iso_;

        //- Per point: originating mesh vertex/cc. See encoding above
        edgeList pointToVerts_;

        //- For every face the original cell in mesh
        labelList meshCells_;

        //- For every point the originating face in mesh
        labelListList pointToFace_;


    // Private Member Functions

        //- Does any edge of triangle cross iso value?
        bool isTriCut
        (
            const triFace& tri,
            const scalarField& pointValues
        ) const;

        //- Determine whether cell is cut
        cellCutType calcCutType
        (
            const bool isTet,
            const label
        ) const;

        //- Determine for all mesh whether cell is cut
        label calcCutTypes
        (
            tetMatcher& tet,
            List<cellCutType>& cellCutTypes
        );

        //- Generate single point on edge
        label generatePoint
        (
            const label facei,
            const bool edgeIsDiag,
            const edge& vertices,

            DynamicList<edge>& pointToVerts,
            DynamicList<labelList>& pointToFace,
            DynamicList<bool>& pointFromDiag,
            EdgeMap<label>& vertsToPoint
        ) const;

        //- Generate triangles from tet
        void generateTriPoints
        (
            const label facei,
            const FixedList<scalar, 4>& s,
            const FixedList<point, 4>& p,
            const FixedList<label, 4>& pIndex,
            const FixedList<bool, 6>& edgeIsDiag,

            DynamicList<edge>& pointToVerts,
            DynamicList<labelList>& pointToFace,
            DynamicList<bool>& pointFromDiag,

            EdgeMap<label>& vertsToPoint,
            DynamicList<label>& verts
        ) const;

        //- Generate triangles from cell
        void generateTriPoints
        (
            const polyMesh& mesh,
            const label celli,
            const bool isTet,

            DynamicList<edge>& pointToVerts,
            DynamicList<labelList>& pointToFace,
            DynamicList<bool>& pointFromDiag,

            EdgeMap<label>& vertsToPoint,
            DynamicList<label>& verts,
            DynamicList<label>& faceLabels
        ) const;


        // Simplification

            label findEdge
            (
                const face& f,
                const edge& e
            ) const;

            label findMeshEdge
            (
                const boolList& pointFromDiag,
                const primitivePatch& pp,
                const labelList& loop,
                const label startFp
            ) const;

            void detectFaceEdges
            (
                const primitivePatch& pp,
                const labelList& loop,
                const boolList& pointFromDiag,
                const labelListList& pointToFace,
                EdgeMap<label>& vertsToFace
            ) const;

            void filterFace
            (
                const primitivePatch& pp,
                const labelList& loop,
                const boolList& pointFromDiag,
                const labelListList& pointToFace,
                face& f,
                const EdgeMap<label>& vertsToFace
            ) const;

            void triangulateOutside
            (
                const bool filterDiag,
                const PrimitivePatch<face, SubList, const pointField&>& pp,
                const boolList& pointFromDiag,
                const labelListList& pointToFace,
                const label cellID,
                const EdgeMap<label>& vertsToFace,
                DynamicList<face>& compactFaces,
                DynamicList<label>& compactCellIDs
            ) const;

            MeshedSurface<face> removeInsidePoints
            (
                const bool filterDiag,
                const MeshedSurface<face>& s,
                const boolList& pointFromDiag,
                const labelListList& pointToFace,
                const labelList& start,              // per cell:starting tri
                const EdgeMap<label>& vertsToFace,
                DynamicList<label>& pointCompactMap, // per point the original
                DynamicList<label>& compactCellIDs   // per face the cellID
            ) const;


public:

    //- Runtime type information
    TypeName("isoSurfaceTopo");


    // Constructors

        //- Construct from dictionary
        isoSurfaceTopo
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const filterType filter = DIAGCELL
        );


    // Member Functions

        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- For every point originating face (pyramid) in mesh
        const labelListList& pointToFace() const
        {
            return pointToFace_;
        }

        //- Per point: originating mesh vertex/cc. See encoding above�
        const edgeList& pointToVerts() const
        {
            return pointToVerts_;
        }

        //- Interpolates cCoords,pCoords.
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const Field<Type>& cCoords,
            const Field<Type>& pCoords
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isoSurfaceTopoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
