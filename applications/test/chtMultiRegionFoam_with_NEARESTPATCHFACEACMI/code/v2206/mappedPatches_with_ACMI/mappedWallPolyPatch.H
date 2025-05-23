/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::mappedWallPolyPatch

Description
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

Note
    Storage is not optimal. It stores all face centres and cells on all
    processors to keep the addressing calculation simple.

SourceFiles
    mappedWallPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedWallPolyPatch_H
#define mappedWallPolyPatch_H

#include "mappedPatchBase.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class mappedWallPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedWallPolyPatch
:
    public wallPolyPatch,
    public mappedPatchBase2
{

protected:

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("mappedWall2");


    // Constructors

        //- Construct from components
        mappedWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components
        mappedWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& sampleRegion,
            const mappedPatchBase2::sampleMode mode,
            const word& samplePatch,
            const vectorField& offset,
            const polyBoundaryMesh& bm
        );

        //- Construct from components. Uniform offset.
        mappedWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& sampleRegion,
            const mappedPatchBase2::sampleMode mode,
            const word& samplePatch,
            const vector& offset,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        mappedWallPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        mappedWallPolyPatch
        (
            const mappedWallPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //- face list and boundary mesh information
        mappedWallPolyPatch
        (
            const mappedWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        mappedWallPolyPatch
        (
            const mappedWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new mappedWallPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //- and boundary mesh
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
                new mappedWallPolyPatch
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
        //- and boundary mesh
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
                new mappedWallPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }

        // Implicit treatment functions

            //- Return number of new internal of this polyPatch faces
            virtual void newInternalProcFaces(label& iFaces, label& pFaces) const
            {
                label nbrFaces = lookupPatch(sampleRegion_, samplePatch_).size();
                iFaces = patch_.size();
                pFaces = max(iFaces - nbrFaces, 0);
            }


            //- Return nbrCells
            virtual const labelUList& nbrCells() const
            {
                return lookupPatch(sampleRegion_, samplePatch_).faceCells();
            }

            //- Return nbr patch ID
            virtual label neighbPolyPatchID() const
            {
                return lookupPatch(sampleRegion_, samplePatch_).index();
            }

            //- Return collocated faces map
            virtual refPtr<labelListList> mapCollocatedFaces() const
            {
                refPtr<labelListList> tMap(new labelListList(patch_.size()));
                labelListList& map = tMap.ref();
                forAll (map, i)
                {
                    labelList& subMap = map[i];
                    subMap.setSize(1);
                    subMap[0] = i;
                }
                return tMap;
            }

            //- Return implicit master
            virtual bool masterImplicit() const
            {
                return owner();
            }

            //- Return neigh region ID
            virtual word neighbRegionID() const
            {
                return sampleRegion_;
            }


    //- Destructor
    virtual ~mappedWallPolyPatch();


    // Member functions

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
