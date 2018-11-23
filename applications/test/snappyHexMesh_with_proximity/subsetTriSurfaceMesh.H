/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::subsetTriSurfaceMesh

Description
    Version of triSurfaceMesh that subsets based on read field

SourceFiles
    subsetTriSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef subsetTriSurfaceMesh_H
#define subsetTriSurfaceMesh_H

#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class subsetTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class subsetTriSurfaceMesh
:
    public triSurfaceMesh
{
    // Private member data

    // Private Member Functions

        static triSurface subsetTriSurface
        (
            const IOobject& io,
            const dictionary& dict
        );


        //- Disallow default bitwise copy construct
        subsetTriSurfaceMesh(const subsetTriSurfaceMesh&);

        //- Disallow default bitwise assignment
        void operator=(const subsetTriSurfaceMesh&);


public:

    //- Runtime type information
    TypeName("subsetTriSurfaceMesh");


    // Constructors

//         //- Construct from triSurface
//         subsetTriSurfaceMesh(const IOobject&, const triSurface&);
// 
//         //- Construct read
//         subsetTriSurfaceMesh(const IOobject& io);
// 
        //- Construct from IO and dictionary (used by searchableSurface).
        //  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
        subsetTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~subsetTriSurfaceMesh();

//         //- Clear storage
//         void clearOut();


    // Member Functions

//         //- Move points
//         virtual void movePoints(const pointField&);
// 
//         //- Demand driven construction of octree for boundary edges
//         const indexedOctree<treeDataEdge>& edgeTree() const;


        // searchableSurface implementation

//             virtual const wordList& regions() const;
// 
//             //- Whether supports volume type below. I.e. whether is closed.
//             virtual bool hasVolumeType() const;
// 
//             //- Range of local indices that can be returned.
//             virtual label size() const
//             {
//                 return triSurface::size();
//             }
// 
//             //- Get representative set of element coordinates
//             //  Usually the element centres (should be of length size()).
//             virtual tmp<pointField> coordinates() const;
// 
//             //- Get bounding spheres (centre and radius squared). Any point
//             //  on surface is guaranteed to be inside.
//             virtual void boundingSpheres
//             (
//                 pointField& centres,
//                 scalarField& radiusSqr
//             ) const;
// 
//             //- Get the points that define the surface.
//             virtual tmp<pointField> points() const;
// 
//             // Does any part of the surface overlap the supplied bound box?
//             virtual bool overlaps(const boundBox& bb) const;
// 
//             virtual void findNearest
//             (
//                 const pointField& sample,
//                 const scalarField& nearestDistSqr,
//                 List<pointIndexHit>&
//             ) const;
// 
//             virtual void findNearest
//             (
//                 const pointField& sample,
//                 const scalarField& nearestDistSqr,
//                 const labelList& regionIndices,
//                 List<pointIndexHit>&
//             ) const;
// 
//             virtual void findLine
//             (
//                 const pointField& start,
//                 const pointField& end,
//                 List<pointIndexHit>&
//             ) const;
// 
//             virtual void findLineAny
//             (
//                 const pointField& start,
//                 const pointField& end,
//                 List<pointIndexHit>&
//             ) const;
// 
//             //- Get all intersections in order from start to end.
//             virtual void findLineAll
//             (
//                 const pointField& start,
//                 const pointField& end,
//                 List<List<pointIndexHit>>&
//             ) const;
// 
//             //- From a set of points and indices get the region
//             virtual void getRegion
//             (
//                 const List<pointIndexHit>&,
//                 labelList& region
//             ) const;
// 
//             //- From a set of points and indices get the normal
//             virtual void getNormal
//             (
//                 const List<pointIndexHit>&,
//                 vectorField& normal
//             ) const;
// 
//             //- Determine type (inside/outside/mixed) for point. unknown if
//             //  cannot be determined (e.g. non-manifold surface)
//             virtual void getVolumeType
//             (
//                 const pointField&,
//                 List<volumeType>&
//             ) const;
// 
// 
//         // Other
// 
//             //- WIP. Store element-wise field.
//             virtual void setField(const labelList& values);
// 
//             //- WIP. From a set of hits (points and
//             //  indices) get the specified field. Misses do not get set.
//             virtual void getField(const List<pointIndexHit>&, labelList&) const;
// 
//             //- Return a pair of triSurfaceScalarFields representing the
//             //  internal and external closeness of regions of the surface
//             Pair<tmp<triSurfaceScalarField>> extractCloseness
//             (
//                 const scalar internalAngleTolerance = 80,
//                 const scalar externalAngleTolerance = 10
//             ) const;
// 
//             //- Return a pair of triSurfaceScalarPointFields representing the
//             //  internal and external closeness of regions of the surface
//             Pair<tmp<triSurfacePointScalarField>> extractPointCloseness
//             (
//                 const scalar internalAngleTolerance = 80,
//                 const scalar externalAngleTolerance = 10
//             ) const;
// 
// 
//         // regIOobject implementation
// 
//             bool writeData(Ostream&) const
//             {
//                 NotImplemented;
//                 return false;
//             }
// 
//             //- Write using given format, version and compression
//             virtual bool writeObject
//             (
//                 IOstream::streamFormat fmt,
//                 IOstream::versionNumber ver,
//                 IOstream::compressionType cmp,
//                 const bool valid
//             ) const;
// 
//             //- Is object global
//             virtual bool global() const
//             {
//                 return true;
//             }
// 
//             //- Return complete path + object name if the file exists
//             //  either in the case/processor or case otherwise null
//             virtual fileName filePath() const
//             {
//                 return searchableSurface::globalFilePath(type());
//             }
};


//- Template function for obtaining global status
template<>
inline bool typeGlobal<subsetTriSurfaceMesh>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
