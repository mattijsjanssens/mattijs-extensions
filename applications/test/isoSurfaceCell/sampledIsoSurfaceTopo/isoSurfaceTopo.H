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

        //- Whether cell might be cut
        List<cellCutType> cellCutType_;

        //- Estimated number of cut cells
        label nCutCells_;

        //- Per point: originating mesh vertex/cc. See encoding above�
        edgeList pointToVerts_;

        //- For every face the original cell in mesh
        labelList meshCells_;

        //- For every point the originating face in mesh
        labelList pointToFace_;


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
        void calcCutTypes(tetMatcher& tet);

        //- Generate single point on edge
        label generatePoint
        (
            const label facei,
            const bool edgeIsDiag,
            const edge& vertices,

            DynamicList<edge>& pointToVerts,
            DynamicList<label>& pointToFace,
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
            DynamicList<label>& pointToFace,
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
            DynamicList<label>& pointToFace,
            DynamicList<bool>& pointFromDiag,

            EdgeMap<label>& vertsToPoint,
            DynamicList<label>& verts,
            DynamicList<label>& faceLabels
        ) const;


        // Simplification

            void triangulateOutside
            (
                const bool filterDiag,
                const PrimitivePatch<face, SubList, const pointField&>& pp,
                const boolList& pointFromDiag,
                const labelList& pointToFace,
                const label cellID,

                DynamicList<face>& compactFaces,
                DynamicList<label>& compactCellIDs
            ) const;

            MeshedSurface<face> removeInsidePoints
            (
                const bool filterDiag,
                const MeshedSurface<face>& s,
                const boolList& pointFromDiag,
                const labelList& pointToFace,
                const labelList& start,              // per cell:starting tri
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
            const bool regularise,
            const bool removeDiagPoints = true
        );


    // Member Functions

        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- For every point originating face (pyramid) in mesh
        const labelList& pointToFace() const
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