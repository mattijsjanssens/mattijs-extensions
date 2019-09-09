/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

Application

Description

Usage

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurfaceMesh.H"
#include "indexedOctree.H"
//#include "treeBoundBox.H"
//#include "PackedBoolList.H"
//#include "unitConversion.H"
//#include "searchableSurfaces.H"
//#include "IOdictionary.H"

// core/api.cpp*
#include "api.h"
#include "parallel.h"
#include "paramset.h"
#include "spectrum.h"
#include "scene.h"
#include "film.h"
#include "medium.h"
#include "stats.h"

// API Additional Headers
#include "accelerators/bvh.h"
#include "accelerators/kdtreeaccel.h"
//#include "cameras/environment.h"
//#include "cameras/orthographic.h"
//#include "cameras/perspective.h"
//#include "cameras/realistic.h"
//#include "filters/box.h"
//#include "filters/gaussian.h"
//#include "filters/mitchell.h"
//#include "filters/sinc.h"
//#include "filters/triangle.h"
//#include "integrators/bdpt.h"
//#include "integrators/directlighting.h"
//#include "integrators/mlt.h"
//#include "integrators/ao.h"
//#include "integrators/path.h"
//#include "integrators/sppm.h"
//#include "integrators/volpath.h"
//#include "integrators/whitted.h"
//#include "lights/diffuse.h"
//#include "lights/distant.h"
//#include "lights/goniometric.h"
//#include "lights/infinite.h"
//#include "lights/point.h"
//#include "lights/projection.h"
//#include "lights/spot.h"
//#include "materials/disney.h"
//#include "materials/fourier.h"
//#include "materials/glass.h"
//#include "materials/hair.h"
//#include "materials/kdsubsurface.h"
//#include "materials/matte.h"
//#include "materials/metal.h"
//#include "materials/mirror.h"
//#include "materials/mixmat.h"
//#include "materials/plastic.h"
//#include "materials/substrate.h"
//#include "materials/subsurface.h"
//#include "materials/translucent.h"
//#include "materials/uber.h"
//#include "samplers/halton.h"
//#include "samplers/maxmin.h"
//#include "samplers/random.h"
//#include "samplers/sobol.h"
//#include "samplers/stratified.h"
//#include "samplers/zerotwosequence.h"
//#include "shapes/cone.h"
//#include "shapes/curve.h"
//#include "shapes/cylinder.h"
//#include "shapes/disk.h"
//#include "shapes/heightfield.h"
//#include "shapes/hyperboloid.h"
//#include "shapes/loopsubdiv.h"
//#include "shapes/nurbs.h"
//#include "shapes/paraboloid.h"
//#include "shapes/sphere.h"
//#include "shapes/triangle.h"
//#include "shapes/plymesh.h"
//#include "textures/bilerp.h"
//#include "textures/checkerboard.h"
//#include "textures/constant.h"
//#include "textures/dots.h"
//#include "textures/fbm.h"
//#include "textures/imagemap.h"
//#include "textures/marble.h"
//#include "textures/mix.h"
//#include "textures/ptex.h"
//#include "textures/scale.h"
//#include "textures/uv.h"
//#include "textures/windy.h"
//#include "textures/wrinkled.h"
//#include "media/grid.h"
//#include "media/homogeneous.h"

#include <map>
#include <stdio.h>



using namespace Foam;

//// Split facei along edgeI at position newPointi
//void greenRefine
//(
//    const triSurface& surf,
//    const label facei,
//    const label edgeI,
//    const label newPointi,
//    DynamicList<labelledTri>& newFaces
//)
//{
//    const labelledTri& f = surf.localFaces()[facei];
//    const edge& e = surf.edges()[edgeI];
//
//    // Find index of edge in face.
//
//    label fp0 = findIndex(f, e[0]);
//    label fp1 = f.fcIndex(fp0);
//    label fp2 = f.fcIndex(fp1);
//
//    if (f[fp1] == e[1])
//    {
//        // Edge oriented like face
//        newFaces.append
//        (
//            labelledTri
//            (
//                f[fp0],
//                newPointi,
//                f[fp2],
//                f.region()
//            )
//        );
//        newFaces.append
//        (
//            labelledTri
//            (
//                newPointi,
//                f[fp1],
//                f[fp2],
//                f.region()
//            )
//        );
//    }
//    else
//    {
//        newFaces.append
//        (
//            labelledTri
//            (
//                f[fp2],
//                newPointi,
//                f[fp1],
//                f.region()
//            )
//        );
//        newFaces.append
//        (
//            labelledTri
//            (
//                newPointi,
//                f[fp0],
//                f[fp1],
//                f.region()
//            )
//        );
//    }
//}
//
//
//void createBoundaryEdgeTrees
//(
//    const PtrList<triSurfaceMesh>& surfs,
//    PtrList<indexedOctree<treeDataEdge>>& bEdgeTrees,
//    labelListList& treeBoundaryEdges
//)
//{
//    forAll(surfs, surfI)
//    {
//        const triSurface& surf = surfs[surfI];
//
//        // Boundary edges
//        treeBoundaryEdges[surfI] =
//            labelList
//            (
//                identity(surf.nEdges() - surf.nInternalEdges())
//              + surf.nInternalEdges()
//            );
//
//        Random rndGen(17301893);
//
//        // Slightly extended bb. Slightly off-centred just so on symmetric
//        // geometry there are less face/edge aligned items.
//        treeBoundBox bb
//        (
//            treeBoundBox(UList<point>(surf.localPoints())).extend(1e-4)
//        );
//
//        bEdgeTrees.set
//        (
//            surfI,
//            new indexedOctree<treeDataEdge>
//            (
//                treeDataEdge
//                (
//                    false,                      // cachebb
//                    surf.edges(),               // edges
//                    surf.localPoints(),         // points
//                    treeBoundaryEdges[surfI]    // selected edges
//                ),
//                bb,     // bb
//                8,      // maxLevel
//                10,     // leafsize
//                3.0     // duplicity
//            )
//       );
//    }
//}
//
//
//class findNearestOpSubset
//{
//    const indexedOctree<treeDataEdge>& tree_;
//
//    DynamicList<label>& shapeMask_;
//
//public:
//
//    findNearestOpSubset
//    (
//        const indexedOctree<treeDataEdge>& tree,
//        DynamicList<label>& shapeMask
//    )
//    :
//        tree_(tree),
//        shapeMask_(shapeMask)
//    {}
//
//    void operator()
//    (
//        const labelUList& indices,
//        const point& sample,
//
//        scalar& nearestDistSqr,
//        label& minIndex,
//        point& nearestPoint
//    ) const
//    {
//        const treeDataEdge& shape = tree_.shapes();
//
//        forAll(indices, i)
//        {
//            const label index = indices[i];
//            const label edgeIndex = shape.edgeLabels()[index];
//
//            if
//            (
//                !shapeMask_.empty()
//             && findIndex(shapeMask_, edgeIndex) != -1
//            )
//            {
//                continue;
//            }
//
//            const edge& e = shape.edges()[edgeIndex];
//
//            pointHit nearHit = e.line(shape.points()).nearestDist(sample);
//
//            // Only register hit if closest point is not an edge point
//            if (nearHit.hit())
//            {
//                scalar distSqr = sqr(nearHit.distance());
//
//                if (distSqr < nearestDistSqr)
//                {
//                    nearestDistSqr = distSqr;
//                    minIndex = index;
//                    nearestPoint = nearHit.rawPoint();
//                }
//            }
//        }
//    }
//};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    const triSurfaceMesh surf
    (
        IOobject
        (
            "junk.obj",
            runTime.constant(),
            "triSurface",
            runTime
        )
    );

    std::shared_ptr<Primitive> accel
    (
        CreateKdTreeAccelerator
        (
            std::move(prims),
            paramSet
        )
    );


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
