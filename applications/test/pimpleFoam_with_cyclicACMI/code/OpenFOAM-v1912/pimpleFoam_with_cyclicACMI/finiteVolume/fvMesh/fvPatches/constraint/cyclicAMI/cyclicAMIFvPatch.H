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
    Foam::cyclicAMIFvPatch

Description
    Cyclic patch for Arbitrary Mesh Interface (AMI)

SourceFiles
    cyclicAMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIFvPatch_H
#define cyclicAMIFvPatch_H

#include "cyclicAMIPolyPatch.H"
#include "coupledFvPatch.H"
#include "cyclicAMILduInterface.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicAMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIFvPatch
:
    public coupledFvPatch,
    public cyclicAMILduInterface
{
    // Private data

        const cyclicAMIPolyPatch& cyclicAMIPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;


public:

    //- Runtime type information
    TypeName(cyclicAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicAMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            cyclicAMILduInterface(),
            cyclicAMIPolyPatch_(refCast<const cyclicAMIPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const cyclicAMIPolyPatch& cyclicAMIPatch() const
            {
                return cyclicAMIPolyPatch_;
            }

            //- Return all neighbours (including ones without AMI)
            virtual const labelList& neighbPatchIDs() const
            {
                return cyclicAMIPolyPatch_.neighbPatchIDs();
            }

            //- Am I owner of this 'connection'
            virtual bool owner(const label index) const
            {
                return cyclicAMIPolyPatch_.owner(index);
            }

            //- Return neighbour
            virtual const cyclicAMIFvPatch& neighbPatch(const label index) const
            {
                return refCast<const cyclicAMIFvPatch>
                (
                    this->boundaryMesh()
                    [
                        cyclicAMIPolyPatch_.neighbPatchIDs()[index]
                    ]
                );
            }

            //- Return a reference to the AMI interpolator (or null)
            virtual tmpNrc<AMIPatchToPatchInterpolation> AMI
            (
                const label index
            ) const
            {
                return cyclicAMIPolyPatch_.AMI(index);
            }

            //- Total neighbour size (of neighbours with valid AMI only)
            virtual label neighbSize() const
            {
                return cyclicAMIPolyPatch_.neighbSize();
            }

            //- Indices of neighbours with valid AMI
            virtual const labelList& AMIIndices() const
            {
                return cyclicAMIPolyPatch_.AMIIndices();
            }

            //- Index on neighbour interface of my interface
            virtual label neighbIndex(const label index) const
            {
                return cyclicAMIPolyPatch_.neighbIndex(index);
            }

            //- Return true if applying the low weight correction
            virtual bool applyLowWeightCorrection() const
            {
                return cyclicAMIPolyPatch_.applyLowWeightCorrection();
            }

            //- Are the cyclic planes parallel
            virtual bool parallel() const
            {
                return cyclicAMIPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicAMIPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicAMIPolyPatch_.reverseT();
            }

            const cyclicAMIFvPatch& neighbFvPatch(const label index) const
            {
                return refCast<const cyclicAMIFvPatch>
                (
                    this->boundaryMesh()
                    [
                        cyclicAMIPolyPatch_.neighbPatchIDs()[index]
                    ]
                );
            }

            //- Return true if this patch is coupled. This is equivalent
            //  to the coupledPolyPatch::coupled() if parallel running or
            //  both sides present, false otherwise
            virtual bool coupled() const;

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return cyclicAMIPolyPatch_.interpolate(fld, defaultValues);
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return cyclicAMIPolyPatch_.interpolate(tFld, defaultValues);
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
