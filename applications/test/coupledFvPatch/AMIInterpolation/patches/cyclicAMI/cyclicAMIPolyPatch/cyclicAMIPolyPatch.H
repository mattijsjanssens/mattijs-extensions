/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::cyclicAMIPolyPatch

Description
    Cyclic patch for Arbitrary Mesh Interface (AMI)

SourceFiles
    cyclicAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIPolyPatch_H
#define cyclicAMIPolyPatch_H

#include "coupledPolyPatch.H"
#include "AMIPatchToPatchInterpolation.H"
#include "polyBoundaryMesh.H"
#include "coupleGroupIdentifier.H"

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
        mutable word nbrPatchName_;

        //- Optional patchGroup to find neighbPatch
        const coupleGroupIdentifier coupleGroup_;

        //- Index of other half
        mutable label nbrPatchID_;


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


        //- AMI interpolation classes
        mutable PtrList<AMIPatchToPatchInterpolation> AMIs_;

        //- AMI transforms (from source to target)
        mutable List<vectorTensorTransform> AMITransforms_;

        //- Flag to indicate that slave patch should be reversed for AMI
        const bool AMIReverse_;

        //- Flag to indicate that patches should match/overlap
        const bool AMIRequireMatch_;

        //- Low weight correction threshold for AMI
        const scalar AMILowWeightCorrection_;

        //- AMI Method
        const AMIPatchToPatchInterpolation::interpolationMethod AMIMethod_;

        //- Projection surface
        mutable autoPtr<searchableSurface> surfPtr_;

        //- Dictionary used during projection surface construction
        const dictionary surfDict_;


    // Protected Member Functions

        //- Reset the AMI interpolator
        virtual void resetAMI() const;

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

        //- Construct from (base couped patch) components
        cyclicAMIPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType transform = UNKNOWN,
            const bool AMIRequireMatch = true,
            const AMIPatchToPatchInterpolation::interpolationMethod AMIMethod =
                AMIPatchToPatchInterpolation::imFaceAreaWeight
        );

        //- Construct from dictionary
        cyclicAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const bool AMIRequireMatch = true,
            const AMIPatchToPatchInterpolation::interpolationMethod AMIMethod =
                AMIPatchToPatchInterpolation::imFaceAreaWeight
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
            const word& nbrPatchName
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
            return autoPtr<polyPatch>(new cyclicAMIPolyPatch(*this, bm));
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
                    nbrPatchName_
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

            //- Neighbour patch name
            inline const word& neighbPatchName() const;

            //- Neighbour patch ID
            virtual label neighbPatchID() const;

            //- Does this side own the patch?
            virtual bool owner() const;

            //- Return a reference to the neighbour patch
            virtual const cyclicAMIPolyPatch& neighbPatch() const;

            //- Return a reference to the projection surface
            const autoPtr<searchableSurface>& surfPtr() const;

            //- Return a reference to the AMI interpolators
            const PtrList<AMIPatchToPatchInterpolation>& AMIs() const;

            //- Return a reference to the AMI transforms
            const List<vectorTensorTransform>& AMITransforms() const;

            //- Return true if applying the low weight correction
            bool applyLowWeightCorrection() const;

            //- Return the weights sum for this patch
            virtual const scalarField& weightsSum() const;

            //- Return the weights sum for the neighbour patch
            virtual const scalarField& neighbWeightsSum() const;



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

                //- Interpolate field component
                tmp<scalarField> interpolate
                (
                    const scalarField& field,
                    const direction cmpt,
                    const direction rank,
                    const scalarUList& defaultValues = scalarUList()
                ) const;

//XXXX
                //- Interpolate field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const Field<Type>& fld,
                    const Field<typename outerProduct<vector, Type>::type>&,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;
//XXXX


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

        //- Return the transform and face indices on neighbour patch which
        //  shares point p following trajectory vector n
        labelPair pointAMIAndFace
        (
            const label facei,
            const vector& n,
            point& p
        ) const;

        //- Index of processor that holds all of both sides, or -1 if
        //  distributed
        label singlePatchProc() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};

template<>
Foam::tmp<Foam::Field<Foam::scalar>> Foam::cyclicAMIPolyPatch::interpolate
(
    const Field<scalar>& fld,
    const Field<typename outerProduct<vector, scalar>::type>& gradFld,
    const UList<scalar>& defaultValues
) const;

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
