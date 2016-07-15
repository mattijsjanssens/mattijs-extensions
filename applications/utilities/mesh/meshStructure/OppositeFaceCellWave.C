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

\*---------------------------------------------------------------------------*/

#include "OppositeFaceCellWave.H"
#include "polyMesh.H"
// #include "processorPolyPatch.H"
// #include "cyclicPolyPatch.H"
// #include "cyclicAMIPolyPatch.H"
// #include "OPstream.H"
// #include "IPstream.H"
// #include "PstreamReduceOps.H"
// #include "debug.H"
// #include "typeInfo.H"
// #include "SubField.H"
// #include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class TrackingData>
void Foam::OppositeFaceCellWave<Type, TrackingData>::opposingFaceLabels
(
    const label celli,
    const label masterFaceLabel,
    DynamicList<label>& oppositeFaceLabels
) const
{
    // Variant of cell::opposingFaceLabel

    // Algorithm:
    // Go through all the faces of the cell and find the one which
    // does not share a single vertex with the master face.  If there
    // are two or more such faces, return the first one and issue a
    // warning; if there is no opposite face, return -1;

    const face& masterFace = this->mesh_.faces()[masterFaceLabel];

    const labelList& curFaceLabels = this->mesh_.cells()[celli];

    oppositeFaceLabels.clear();

    forAll(curFaceLabels, facei)
    {
        // Compare the face with the master
        const face& curFace = this->mesh_.faces()[curFaceLabels[facei]];

        // Skip the master face
        if (curFaceLabels[facei] != masterFaceLabel)
        {
            bool sharedPoint = false;

            // Compare every vertex of the current face against the
            // vertices of the master face
            forAll(curFace, pointi)
            {
                const label l = curFace[pointi];

                forAll(masterFace, masterPointi)
                {
                    if (masterFace[masterPointi] == l)
                    {
                        sharedPoint = true;
                        break;
                    }
                }

                if (sharedPoint) break;
            }

            // If no points are shared, this is the opposite face
            if (!sharedPoint)
            {
                // Found opposite face
                oppositeFaceLabels.append(curFaceLabels[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Iterate, propagating changedFacesInfo across mesh, until no change (or
// maxIter reached). Initial cell values specified.
template<class Type, class TrackingData>
Foam::OppositeFaceCellWave<Type, TrackingData>::OppositeFaceCellWave
(
    const polyMesh& mesh,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    const label maxIter,
    TrackingData& td
)
:
    FaceCellWave<Type, TrackingData>
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        allFaceInfo,
        allCellInfo,
        0,              //maxIter,
        td
    ),
    changedOppositeFaces_(this->mesh_.nCells())
{
    // Iterate until nothing changes
    label iter = this->iterate(maxIter);

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter."
            << endl
            << "    maxIter:" << maxIter << endl
            << "    nChangedCells:" << this->changedCells_.size() << endl
            << "    nChangedFaces:" << this->changedFaces_.size() << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type, class TrackingData>
// Foam::labelPair
// Foam::OppositeFaceCellWave<Type, TrackingData>::faceToCellToFace()
// {
//     const labelList& owner = this->mesh_.faceOwner();
//     const labelList& neighbour = this->mesh_.faceNeighbour();
//     label nInternalFaces = this->mesh_.nInternalFaces();
// 
//     DynamicList<label> oppositeFaceLabels;
// 
//     DynamicList<label> newChangedFaces(this->changedFaces_.capacity());
//     PackedBoolList newChangedFace(this->mesh_.nFaces());
// 
//     this->changedCells_.clear();
//     this->changedCell_ = false;
// 
//     forAll(this->changedFaces_, changedFacei)
//     {
//         label facei = this->changedFaces_[changedFacei];
//         Pout<< "face:" << facei << " at:" << this->mesh_.faceCentres()[facei]
//             << " distance:" << this->allFaceInfo_[facei].distance() << endl;
// 
//         if (!this->changedFace_[facei])
//         {
//             FatalErrorInFunction
//                 << "Face " << facei
//                 << " not marked as having been changed"
//                 << abort(FatalError);
//         }
// 
// 
//         const Type& neighbourWallInfo = this->allFaceInfo_[facei];
// 
//         // Evaluate all connected cells
// 
//         // Owner
//         {
//             label celli = owner[facei];
//             Type& currentWallInfo = this->allCellInfo_[celli];
// 
//             if (!currentWallInfo.equal(neighbourWallInfo, this->td_))
//             {
//                 opposingFaceLabels(celli, facei, oppositeFaceLabels);
//                 Pout<< "    celli:" << celli
//                     << " at:" << this->mesh_.cellCentres()[celli]
//                     << " oppositefaces:"
//                     << pointField(this->mesh_.faceCentres(), oppositeFaceLabels)
//                     << endl;
// 
//                 if (oppositeFaceLabels.size())
//                 {
//                     bool propagate = this->updateCell
//                     (
//                         celli,
//                         facei,
//                         neighbourWallInfo,
//                         this->propagationTol_,
//                         currentWallInfo
//                     );
// 
//                     if (propagate && oppositeFaceLabels.size() == 1)
//                     {
//                         const Type& neighbourWallInfo =
//                             this->allCellInfo_[celli];
// 
//                         label oppFacei = oppositeFaceLabels[0];
// 
//                         Type& currentWallInfo = this->allFaceInfo_[oppFacei];
// 
//                         if
//                         (
//                            !currentWallInfo.equal(neighbourWallInfo, this->td_)
//                         )
//                         {
//                             bool propagate = this->updateFace
//                             (
//                                 oppFacei,
//                                 celli,
//                                 neighbourWallInfo,
//                                 this->propagationTol_,
//                                 currentWallInfo
//                             );
//                             if (propagate)
//                             {
//                                 newChangedFaces.append(oppFacei);
//                                 newChangedFace[oppFacei] = true;
//                             }
//                         }
//                     }
//                 }
//             }
//         }
// 
//         // Neighbour
//         if (facei < nInternalFaces)
//         {
//             label celli = neighbour[facei];
//             Type& currentWallInfo2 = this->allCellInfo_[celli];
// 
//             if (!currentWallInfo2.equal(neighbourWallInfo, this->td_))
//             {
//                 opposingFaceLabels(celli, facei, oppositeFaceLabels);
//                 Pout<< "    nei:" << celli
//                     << " at:" << this->mesh_.cellCentres()[celli]
//                     << " oppositefaces:"
//                     << pointField(this->mesh_.faceCentres(), oppositeFaceLabels)
//                     << endl;
// 
//                 if (oppositeFaceLabels.size())
//                 {
//                     bool propagate = this->updateCell
//                     (
//                         celli,
//                         facei,
//                         neighbourWallInfo,
//                         this->propagationTol_,
//                         currentWallInfo2
//                     );
// 
//                     if (propagate && oppositeFaceLabels.size() == 1)
//                     {
//                         const Type& neighbourWallInfo =
//                             this->allCellInfo_[celli];
// 
//                         label oppFacei = oppositeFaceLabels[0];
// 
//                         Type& currentWallInfo = this->allFaceInfo_[oppFacei];
// 
//                         if
//                         (
//                            !currentWallInfo.equal(neighbourWallInfo, this->td_)
//                         )
//                         {
//                             bool propagate = this->updateFace
//                             (
//                                 oppFacei,
//                                 celli,
//                                 neighbourWallInfo,
//                                 this->propagationTol_,
//                                 currentWallInfo
//                             );
//                             if (propagate)
//                             {
//                                 newChangedFaces.append(oppFacei);
//                                 newChangedFace[oppFacei] = true;
//                             }
//                         }
//                     }
//                 }
//             }
//         }
// 
//         // Reset status of face
//         this->changedFace_[facei] = false;
//     }
// 
//     // Update changed faces
//     this->changedFaces_.transfer(newChangedFaces);
//     this->changedFace_ = newChangedFace;
// 
// 
//     //if (debug & 2)
//     {
//         Pout<< " Changed faces            : " << this->changedFaces_.size()
//             << endl;
//         Pout<< " Changed cells            : " << this->changedCells_.size()
//             << endl;
//     }
// 
//     labelPair changed
//     (
//         returnReduce(this->changedCells_.size(), sumOp<label>()),
//         returnReduce(this->changedFaces_.size(), sumOp<label>())
//     );
// 
//     if (this->hasCyclicPatches_)
//     {
//         // Transfer changed faces across cyclic halves
//         this->handleCyclicPatches();
//     }
// 
//     if (this->hasCyclicAMIPatches_)
//     {
//         this->handleAMICyclicPatches();
//     }
// 
//     if (Pstream::parRun())
//     {
//         // Transfer changed faces from neighbouring processors.
//         this->handleProcPatches();
//     }
// 
//     return changed;
// }


template<class Type, class TrackingData>
Foam::label Foam::OppositeFaceCellWave<Type, TrackingData>::faceToCell()
{
    const labelList& owner = this->mesh_.faceOwner();
    const labelList& neighbour = this->mesh_.faceNeighbour();
    label nInternalFaces = this->mesh_.nInternalFaces();

    DynamicList<label> oppositeFaceLabels;

    forAll(this->changedFaces_, changedFacei)
    {
        label facei = this->changedFaces_[changedFacei];

        Pout<< "face:" << facei << " at:" << this->mesh_.faceCentres()[facei]
            << " distance:" << this->allFaceInfo_[facei].distance() << endl;


        if (!this->changedFace_[facei])
        {
            FatalErrorInFunction
                << "Face " << facei
                << " not marked as having been changed"
                << abort(FatalError);
        }


        const Type& neighbourWallInfo = this->allFaceInfo_[facei];

        // Evaluate all connected cells

        // Owner
        {
            label celli = owner[facei];
            Type& currentWallInfo = this->allCellInfo_[celli];

            if (!currentWallInfo.equal(neighbourWallInfo, this->td_))
            {
                // Check if cell is prismatic w.r.t facei
                opposingFaceLabels(celli, facei, oppositeFaceLabels);

                Pout<< "    celli:" << celli
                    << " at:" << this->mesh_.cellCentres()[celli]
                    << " oppositefaces:"
                    << pointField(this->mesh_.faceCentres(), oppositeFaceLabels)
                    << endl;

                if (oppositeFaceLabels.size())
                {
                    label sz = this->changedCells_.size();
                    this->updateCell
                    (
                        celli,
                        facei,
                        neighbourWallInfo,
                        this->propagationTol_,
                        currentWallInfo
                    );
                    if
                    (
                        this->changedCells_.size() > sz
                     && oppositeFaceLabels.size() == 1
                    )
                    {
                        changedOppositeFaces_.append(oppositeFaceLabels[0]);
                    }
                }
            }
        }

        // Neighbour.
        if (facei < nInternalFaces)
        {
            label celli = neighbour[facei];
            Type& currentWallInfo2 = this->allCellInfo_[celli];

            if (!currentWallInfo2.equal(neighbourWallInfo, this->td_))
            {
                // Check if cell is prismatic w.r.t facei
                opposingFaceLabels(celli, facei, oppositeFaceLabels);

                Pout<< "    celli:" << celli
                    << " at:" << this->mesh_.cellCentres()[celli]
                    << " oppositefaces:"
                    << pointField(this->mesh_.faceCentres(), oppositeFaceLabels)
                    << endl;

                if (oppositeFaceLabels.size())
                {
                    label sz = this->changedCells_.size();
                    this->updateCell
                    (
                        celli,
                        facei,
                        neighbourWallInfo,
                        this->propagationTol_,
                        currentWallInfo2
                    );
                    if
                    (
                        this->changedCells_.size() > sz
                     && oppositeFaceLabels.size() == 1
                    )
                    {
                        changedOppositeFaces_.append(oppositeFaceLabels[0]);
                    }
                }
            }
        }

        // Reset status of face
        this->changedFace_[facei] = false;
    }

    // Handled all changed faces by now
    this->changedFaces_.clear();

    //if (debug & 2)
    {
        Pout<< " Changed cells            : " << this->changedCells_.size()
            << endl;
    }

    // Sum changedCells over all procs
    label totNChanged = this->changedCells_.size();

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


template<class Type, class TrackingData>
Foam::label Foam::OppositeFaceCellWave<Type, TrackingData>::cellToFace()
{
    forAll(this->changedCells_, changedCelli)
    {
        label celli = this->changedCells_[changedCelli];
        label facei = changedOppositeFaces_[changedCelli];

        Pout<< "cell:" << celli << " at:" << this->mesh_.cellCentres()[celli]
            << " distance:" << this->allCellInfo_[celli].distance()
            << " with oppface:" << facei
            << " at:" << this->mesh_.faceCentres()[facei]
            << endl;


        if (!this->changedCell_[celli])
        {
            FatalErrorInFunction
                << "Cell " << celli << " not marked as having been changed"
                << abort(FatalError);
        }

        const Type& neighbourWallInfo = this->allCellInfo_[celli];

        // Evaluate facei

        Type& currentWallInfo = this->allFaceInfo_[facei];

        if (!currentWallInfo.equal(neighbourWallInfo, this->td_))
        {
            this->updateFace
            (
                facei,
                celli,
                neighbourWallInfo,
                this->propagationTol_,
                currentWallInfo
            );
        }

        // Reset status of cell
        this->changedCell_[celli] = false;
    }

    // Handled all changed cells by now
    this->changedCells_.clear();
    changedOppositeFaces_.clear();

    if (this->hasCyclicPatches_)
    {
        // Transfer changed faces across cyclic halves
        this->handleCyclicPatches();
    }

    if (this->hasCyclicAMIPatches_)
    {
        this->handleAMICyclicPatches();
    }

    if (Pstream::parRun())
    {
        // Transfer changed faces from neighbouring processors.
        this->handleProcPatches();
    }

    //if (debug & 2)
    {
        Pout<< " Changed faces            : " << this->changedFaces_.size()
            << endl;
    }

    // Sum nChangedFaces over all procs
    label totNChanged = this->changedFaces_.size();

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// template<class Type, class TrackingData>
// Foam::label
// Foam::OppositeFaceCellWave<Type, TrackingData>::iterate(const label maxIter)
// {
// DebugVar(maxIter);
// 
//     if (this->hasCyclicPatches_)
//     {
//         // Transfer changed faces across cyclic halves
//         this->handleCyclicPatches();
//     }
// 
//     if (this->hasCyclicAMIPatches_)
//     {
//         this->handleAMICyclicPatches();
//     }
// 
//     if (Pstream::parRun())
//     {
//         // Transfer changed faces from neighbouring processors.
//         this->handleProcPatches();
//     }
// 
//     label iter = 0;
// 
//     while (iter < maxIter)
//     {
//         if (debug)
//         {
//             Info<< " Iteration " << iter << endl;
//         }
// 
//         this->nEvals_ = 0;
// 
//         labelPair changed = faceToCellToFace();
// 
//         //if (debug)
//         {
//             Info<< " Total changed cells      : " << changed.first() << nl
//                 << " Total changed faces      : " << changed.second() << nl
//                 << " Total evaluations        : " << this->nEvals_ << nl
//                 << " Remaining unvisited cells: " << this->nUnvisitedCells_
//                 << nl
//                 << " Remaining unvisited faces: " << this->nUnvisitedFaces_
//                 << endl;
//         }
// 
//         if (changed.second() == 0)
//         {
//             break;
//         }
// 
//         ++iter;
//     }
// 
//     return iter;
// }


// ************************************************************************* //
