/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::cyclicAMIPolyPatch

Description
    Cyclic patch for Arbitrary Mesh Interface (AMI)

    E.g. in polyMesh/boundary
    leftAMI
    {
        type            cyclicAMI;
        inGroups        (coupleGroup);
        nFaces          4;
        startFace       40;
        coupleGroup     coupleGroup;
        method          partialFaceAreaWeightAMI;
    }

    - the coupleGroup (or neighbourPatches) specifies which other patches it
      might talk to
    - if so there is an AMI on either *this or the neighbourPatch, depending
      on which is the lower numbered one.
    - the coupleGroup might have lots and lots of entries but only some
      are filled.

SourceFiles
    cyclicAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIPolyPatch_H
#define cyclicAMIPolyPatch_H

#include "coupleGroupIdentifier.H"
#include "coupledPolyPatch.H"
#include "AMIPatchToPatchInterpolation.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicAMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIPolyPatch
:
    public coupledPolyPatch
{
    // Private Member Functions

        //- Return normal of face at max distance from rotation axis
        vector findFaceNormalMaxRadius(const pointField& faceCentres) const;

        void calcTransforms
        (
            const primitivePatch& half0,
            const pointField& half0Ctrs,
            const vectorField& half0Areas,
            const pointField& half1Ctrs,
            const vectorField& half1Areas
        );


protected:

    // Protected data

        //- Name of other half
        mutable wordList nbrPatchNames_;

        //- Optional patchGroup to find neighbPatch
        const coupleGroupIdentifier coupleGroup_;

        //- Index of other half
        mutable labelList nbrPatchIDs_;


        // Transformations

            // For rotation

                //- Axis of rotation for rotational cyclics
                vector rotationAxis_;

                //- Point on axis of rotation for rotational cyclics
                point rotationCentre_;

                //- Flag to show whether the rotation angle is defined
                bool rotationAngleDefined_;

                //- Rotation angle
                scalar rotationAngle_;


            // For translation

                //- Translation vector
                vector separationVector_;


        //- AMI interpolation class
        mutable PtrList<AMIPatchToPatchInterpolation> AMIPtrs_;

        //- AMI method
        const AMIPatchToPatchInterpolation::interpolationMethod AMIMethod_;

        //- Flag to indicate that slave patch should be reversed for AMI
        const bool AMIReverse_;

        //- Flag to indicate that patches should match/overlap
        bool AMIRequireMatch_;

        //- Low weight correction threshold for AMI
        const scalar AMILowWeightCorrection_;

        //- Indices of neighbours with a valid AMI
        mutable labelList AMIIndices_;

        //- Pre-calculated index of AMI in neighbour
        mutable labelList neighbIndex_;

        //- Overall weight for all *this faces (only when using multiple AMI;
        //  unsized if not)
        mutable scalarField weightsSum_;

        //- Overall size
        mutable label neighbSize_;

        //- Projection surface
        mutable autoPtr<searchableSurface> surfPtr_;

        //- Dictionary used during projection surface construction
        const dictionary surfDict_;


    // Protected Member Functions

        //- Reset the AMI interpolator
        virtual void resetAMI
        (
            const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod =
                AMIPatchToPatchInterpolation::imFaceAreaWeight
        ) const;

        //- Recalculate the transformation tensors
        virtual void calcTransforms();

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


public:

    //- Runtime type information
    TypeName("cyclicAMI");


    // Constructors

        //- Construct from (base coupled patch) components
        cyclicAMIPolyPatch
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
        cyclicAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicAMIPolyPatch(const cyclicAMIPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicAMIPolyPatch
        (
            const cyclicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const wordList& nbrPatchNames
        );

        //- Construct given the original patch and a map
        cyclicAMIPolyPatch
        (
            const cyclicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>::NewFrom<cyclicAMIPolyPatch>(*this, bm);
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
                new cyclicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchNames_
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
                new cyclicAMIPolyPatch
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
    virtual ~cyclicAMIPolyPatch();


    // Member Functions

        // Access

            //- Is patch 'coupled'. Note that on AMI the geometry is not
            //  coupled but the fields are!
            virtual bool coupled() const
            {
                return false;
            }

            //- Neighbour patch names
            inline const wordList& neighbPatchNames() const;

            //- Neighbour patch ID
            virtual const labelList& neighbPatchIDs() const;

            //- Indices of neighbours with valid AMI
            virtual const labelList& AMIIndices() const;

            //- Index on neighbour of AMI
            virtual label neighbIndex(const label index) const;

            //- Does this side own the patch compared to neighbour index?
            virtual bool owner(const label index) const;

            //- Total neighbour size (of neighbours with valid AMI)
            virtual label neighbSize() const;

            //- Return neighbour internal cell data (in AMI/neighbPatchIDs order
            //  compatible with interpolate routines)
            template<class Type>
            tmp<Field<Type>> patchNeighbourField(const UList<Type>& iF) const;


            //- Backwards compatibility

                //- Neighbour patch name
                inline const word& neighbPatchName() const;

                //- Neighbour patch ID
                virtual label neighbPatchID() const;

                //- Return a reference to the neighbour patch
                virtual const cyclicAMIPolyPatch& neighbPatch
                (
                    const label index
                ) const;

                //- Does this side own the patch (only for single neighbour)
                virtual bool owner() const
                {
                    return owner(0);
                }


            //- Return a reference to the projection surface
            const autoPtr<searchableSurface>& surfPtr() const;

            //- Return a reference to the AMI interpolator
            tmpNrc<AMIPatchToPatchInterpolation> AMI(const label index) const;

            //- Return true if applying the low weight correction
            bool applyLowWeightCorrection() const;

            //- Return sum of weights on *this
            const scalarField& weightsSum() const;


            // Transformations

                //- Axis of rotation for rotational cyclic AMI
                inline const vector& rotationAxis() const;

                //- Point on axis of rotation for rotational cyclic AMI
                inline const point& rotationCentre() const;

                //- Translation vector for translational cyclic AMI
                inline const vector& separationVector() const;

                //- Transform patch-based positions from nbr side to this side
                virtual void transformPosition(pointField&) const;

                //- Transform a patch-based position from nbr side to this side
                virtual void transformPosition
                (
                    point& l,
                    const label facei
                ) const;

                //- Transform a patch-based position from this side to nbr side
                virtual void reverseTransformPosition
                (
                    point& l,
                    const label facei
                ) const;

                //- Transform a patch-based direction from this side to nbr side
                virtual void reverseTransformDirection
                (
                    vector& d,
                    const label facei
                ) const;


            // Interpolations

                //- Interpolate field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const Field<Type>& fld,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Interpolate tmp field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const tmp<Field<Type>>& tFld,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;


        //- Calculate the patch geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const pointField& thisCtrs,
            const vectorField& thisAreas,
            const pointField& thisCc,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const pointField& nbrCc
        );

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

        //- Return face index on neighbour patch which shares point p
        //  following trajectory vector n
        label pointFace
        (
            const label facei,
            const vector& n,
            point& p
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cyclicAMIPolyPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicAMIPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
