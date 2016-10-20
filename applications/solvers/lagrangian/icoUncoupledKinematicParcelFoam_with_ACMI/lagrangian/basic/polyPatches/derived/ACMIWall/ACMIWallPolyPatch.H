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
    Foam::ACMIWallPolyPatch

Description
    Wall for use with ACMI to enable Lagrangian support.

SourceFiles
    mappedWallPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef ACMIWallPolyPatch_H
#define ACMIWallPolyPatch_H

#include "wallPolyPatch.H"
#include "coupleGroupIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class cyclicACMIPolyPatch;

/*---------------------------------------------------------------------------*\
                      Class ACMIWallPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class ACMIWallPolyPatch
:
    public wallPolyPatch
{

protected:

        //- PatchGroup to find ACMI patch
        const coupleGroupIdentifier coupleGroup_;

        //- ACMI patch 
        mutable word ACMIPatch_;

        //- ACMI patch index
        mutable label ACMIPatchID_;


//        //- Initialise the calculation of the patch geometry
//        virtual void initGeometry(PstreamBuffers&);
//
//        //- Calculate the patch geometry
//        virtual void calcGeometry(PstreamBuffers&);
//
//        //- Initialise the patches for moving points
//        virtual void initMovePoints(PstreamBuffers&, const pointField&);
//
//        //- Correct patches after moving points
//        virtual void movePoints(PstreamBuffers&, const pointField&);
//
//        //- Initialise the update of the patch topology
//        virtual void initUpdateMesh(PstreamBuffers&);
//
//        //- Update of the patch topology
//        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("ACMIWall");


    // Constructors

        //- Construct from components
        ACMIWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        ACMIWallPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        ACMIWallPolyPatch
        (
            const ACMIWallPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        ACMIWallPolyPatch
        (
            const ACMIWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        ACMIWallPolyPatch
        (
            const ACMIWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new ACMIWallPolyPatch(*this, bm));
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
                new ACMIWallPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
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
                new ACMIWallPolyPatch
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
    virtual ~ACMIWallPolyPatch();


    // Member functions

        //- Access the ACMI partner patch
        const cyclicACMIPolyPatch& ACMI() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
