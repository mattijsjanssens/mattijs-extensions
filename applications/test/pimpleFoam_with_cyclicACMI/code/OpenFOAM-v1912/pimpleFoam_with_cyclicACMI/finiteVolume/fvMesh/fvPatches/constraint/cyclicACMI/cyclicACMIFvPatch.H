/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::cyclicACMIFvPatch

Description
    Cyclic patch for Arbitrarily Coupled Mesh Interface (ACMI)

SourceFiles
    cyclicACMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIFvPatch_H
#define cyclicACMIFvPatch_H

#include "coupledFvPatch.H"
#include "cyclicACMILduInterface.H"
#include "cyclicACMIPolyPatch.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicACMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIFvPatch
:
    public coupledFvPatch,
    public cyclicACMILduInterface
{
    // Private data

        const cyclicACMIPolyPatch& cyclicACMIPolyPatch_;


protected:

    // Protected Member functions

        //- Update the patch areas after AMI update
        void updateAreas() const;

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;


public:

    //- Runtime type information
    TypeName(cyclicACMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicACMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            cyclicACMILduInterface(),
            cyclicACMIPolyPatch_(refCast<const cyclicACMIPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const cyclicACMIPolyPatch& cyclicACMIPatch() const
            {
                return cyclicACMIPolyPatch_;
            }

            //- Return neighbour
            virtual const labelList& neighbPatchIDs() const
            {
                return cyclicACMIPolyPatch_.neighbPatchIDs();
            }

            //- Total neighbour size (of neighbours with valid AMI only)
            virtual label neighbSize() const
            {
                return cyclicACMIPolyPatch_.neighbSize();
            }

            //- Indices of neighbours with valid AMI
            virtual const labelList& AMIIndices() const
            {
                return cyclicACMIPolyPatch_.AMIIndices();
            }

            //- Index on neighbour of AMI
            virtual label neighbIndex(const label index) const
            {
                return cyclicACMIPolyPatch_.neighbIndex(index);
            }

            //- Return neighbour
            virtual label neighbPatchID(const label index) const
            {
                return cyclicACMIPolyPatch_.neighbPatchIDs()[index];
            }

            virtual bool owner(const label index) const
            {
                return cyclicACMIPolyPatch_.owner(index);
            }

            //- Return neighbour fvPatch
            virtual const cyclicACMIFvPatch& neighbPatch
            (
                const label index
            ) const
            {
                return refCast<const cyclicACMIFvPatch>
                (
                    this->boundaryMesh()
                    [
                        cyclicACMIPolyPatch_.neighbPatchIDs()[index]
                    ]
                );
            }

            //- Return neighbour
            virtual label nonOverlapPatchID() const
            {
                return cyclicACMIPolyPatch_.nonOverlapPatchID();
            }

            //- Return non-overlapping fvPatch
            virtual const fvPatch& nonOverlapPatch() const
            {
                return this->boundaryMesh()[nonOverlapPatchID()];
            }

            //- Return a reference to the AMI interpolator
            virtual tmpNrc<AMIPatchToPatchInterpolation> AMI
            (
                const label index
            ) const
            {
                updateAreas();

                return cyclicACMIPolyPatch_.AMI(index);
            }

            //- Are the cyclic planes parallel
            virtual bool parallel() const
            {
                return cyclicACMIPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicACMIPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicACMIPolyPatch_.reverseT();
            }

            const cyclicACMIFvPatch& neighbFvPatch(const label index) const
            {
                return refCast<const cyclicACMIFvPatch>
                (
                    this->boundaryMesh()
                    [
                        cyclicACMIPolyPatch_.neighbPatchIDs()[index]
                    ]
                );
            }

            //- Return true if this patch is coupled. This is equivalent
            //  to the coupledPolyPatch::coupled() if parallel running or
            //  both sides present, false otherwise
            virtual bool coupled() const;

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            //- Interpolate (make sure to have uptodate areas)
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld
            ) const
            {
                updateAreas();

                return
                    cyclicACMIPolyPatch_.cyclicAMIPolyPatch::interpolate
                    (
                        fld
                    );
            }

            //- Interpolate (make sure to have uptodate areas)
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tfld
            ) const
            {
                return interpolate(tfld());
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
