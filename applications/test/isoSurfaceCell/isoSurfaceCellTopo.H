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
    Foam::isoSurfaceCellTopo

Description

SourceFiles
    isoSurfaceCellTopo.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceCellTopo_H
#define isoSurfaceCellTopo_H

#include "labelPair.H"
#include "pointIndexHit.H"
#include "PackedBoolList.H"
#include "MeshedSurface.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceCellTopo Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceCellTopo
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

        //- iso  value
        const scalar iso_;

        //- Whether cell might be cut
        List<cellCutType> cellCutType_;

        //- Estimated number of cut cells
        label nCutCells_;

        //- Per point: originating mesh vertex/cc. See encoding above·
        edgeList pointToVerts_;

        //- For every face the original cell in mesh
        labelList meshCells_;


    // Private Member Functions

        //- Get location of iso value as fraction inbetween s0,s1
        scalar isoFraction
        (
            const scalar s0,
            const scalar s1
        ) const;

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


public:

    //- Runtime type information
    TypeName("isoSurfaceCellTopo");


    // Constructors

        //- Construct from dictionary
        isoSurfaceCellTopo
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const bool regularise,
            const scalar mergeTol = 1e-10    // fraction of bounding box
        );


    // Member Functions

        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
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
    #include "isoSurfaceCellTopoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
