/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::OppositeFaceCellWave

Description
    Wave propagation of information through grid. Every iteration
    information goes through one layer of cells. Templated on information
    that is transferred.

    Handles parallel and cyclics and non-parallel cyclics.

    Note: whether to propagate depends on the return value of Type::update
    which returns true (i.e. propagate) if the value changes by more than a
    certain tolerance.
    This tolerance can be very strict for normal face-cell and parallel
    cyclics (we use a value of 0.01 just to limit propagation of small changes)
    but for non-parallel cyclics this tolerance can be critical and if chosen
    too small can lead to non-convergence.

SourceFiles
    OppositeFaceCellWave.C

\*---------------------------------------------------------------------------*/

#ifndef OppositeFaceCellWave_H
#define OppositeFaceCellWave_H

#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
//class polyMesh;
//class polyPatch;

/*---------------------------------------------------------------------------*\
                        Class OppositeFaceCellWaveName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(OppositeFaceCellWave);


/*---------------------------------------------------------------------------*\
                           Class OppositeFaceCellWave Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class TrackingData = int>
class OppositeFaceCellWave
:
    public FaceCellWave<Type, TrackingData>,
    public OppositeFaceCellWaveName
{
private:

   // Private Member Functions

        void setFaceInfo(const labelList&, const List<Type>&);

//         label faceToCell();
// 
//         label cellToFace();


protected:

        //- For every entry in changedCells (i.e. the cell front) gives
        //  the face that it needs to transfer to
        DynamicList<label> changedOppositeFaces_;


   // Protected Member Functions

        //- Determine 'opposite' faces (= faces not sharing a vertex) on cell
        void opposingFaceLabels
        (
            const label celli,
            const label facei,
            DynamicList<label>&
        ) const;

//         //- Do walk from face to cell and back to face
//         labelPair faceToCellToFace();


public:

    // Constructors

        //- Construct from mesh and list of changed faces with the Type
        //  for these faces. Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        OppositeFaceCellWave
        (
            const polyMesh&,
            const labelList& initialChangedFaces,
            const List<Type>& changedFacesInfo,
            UList<Type>& allFaceInfo,
            UList<Type>& allCellInfo,
            const label maxIter,
            TrackingData& td = FaceCellWave<Type, TrackingData>::dummyTrackData_
        );


    //- Destructor
    virtual ~OppositeFaceCellWave()
    {};


    // Member Functions

            //- Propagate from face to cell. Returns total number of cells
            //  (over all processors) changed.
            virtual label faceToCell();

            //- Propagate from cell to face. Returns total number of faces
            //  (over all processors) changed. (Faces on processorpatches are
            //  counted double)
            virtual label cellToFace();

//             //- Iterate until no changes or maxIter reached.  Returns actual
//             //  number of iterations.
//             virtual label iterate(const label maxIter);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "OppositeFaceCellWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
