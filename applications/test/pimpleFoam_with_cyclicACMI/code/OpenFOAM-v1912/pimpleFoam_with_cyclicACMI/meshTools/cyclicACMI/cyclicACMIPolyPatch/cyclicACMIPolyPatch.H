/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::cyclicACMIPolyPatch

Description
    Cyclic patch for Arbitrarily Coupled Mesh Interface (ACMI)

SourceFiles
    cyclicACMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIPolyPatch_H
#define cyclicACMIPolyPatch_H

#include "cyclicAMIPolyPatch.H"
#include "AMIPatchToPatchInterpolation.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cyclicACMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIPolyPatch
:
    public cyclicAMIPolyPatch
{

private:

    // Private data

        //- Fraction of face area below which face is considered disconnected
        static const scalar tolerance_;

        //- Name of non-overlapping patch
        const word nonOverlapPatchName_;

        //- Index of non-overlapping patch
        mutable label nonOverlapPatchID_;

//        //- Mask/weighting for source patch
//        mutable List<scalarField> srcMasks_;
//
//        //- Mask/weighting for target patch
//        mutable List<scalarField> tgtMasks_;

        //- Precalculated mask from overall weights for all *this faces
        mutable scalarField maskSum_;

        //- Flag to indicate that AMI has been updated
        mutable bool updated_;


protected:

    // Protected Member Functions

        //- Calculate and print some stats
        void writeStats
        (
            const scalarField& wghtsSum,
            const word& origin,
            const word& patchName
        ) const;

        //- Reset the AMI interpolator
        virtual void resetAMI
        (
            const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod =
                AMIPatchToPatchInterpolation::imFaceAreaWeight
        ) const;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers& pBufs, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers& pBufs, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);

        //- Clear geometry
        virtual void clearGeom();

//        //- Return the mask/weighting for the source patch
//        virtual const scalarField& srcMask(const label index) const;
//
//        //- Return the mask/weighting for the target patch
//        virtual const scalarField& tgtMask(const label index) const;


public:

    //- Runtime type information
    TypeName("cyclicACMI");


    // Constructors

        //- Construct from (base coupled patch) components
        cyclicACMIPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType transform = UNKNOWN
        );

        //- Construct from dictionary
        cyclicACMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicACMIPolyPatch
        (
            const cyclicACMIPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicACMIPolyPatch
        (
            const cyclicACMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const wordList& nbrPatchNames,
            const word& nonOverlapPatchName
        );

        //- Construct given the original patch and a map
        cyclicACMIPolyPatch
        (
            const cyclicACMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>::NewFrom<cyclicACMIPolyPatch>(*this, bm);
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicACMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchNames_,
                    nonOverlapPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicACMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~cyclicACMIPolyPatch();


    // Member Functions

        // Access

            //- Reset the updated flag
            inline void setUpdated(bool flag) const;

            //- Return access to the updated flag
            inline bool updated() const;

            //- Return a reference to the neighbour patch
            virtual const cyclicACMIPolyPatch& neighbPatch
            (
                const label index = 0
            ) const;

            //- Non-overlapping patch name
            inline const word& nonOverlapPatchName() const;

            //- Non-overlapping patch ID
            virtual label nonOverlapPatchID() const;

            //- Return a const reference to the non-overlapping patch
            inline const polyPatch& nonOverlapPatch() const;

            //- Return a reference to the non-overlapping patch
            inline polyPatch& nonOverlapPatch();

            //- Mask field where 1 = overlap, 0 = no-overlap
            inline const scalarField& mask() const;

            //- Overlap tolerance
            inline static scalar tolerance();


        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cyclicACMIPolyPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
