/*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | Copyright (C) 2026 Mattijs Janssens
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
    cgalMeshToPolyMesh

Description
    Test application that reads geometry using CGAL, extracts features,
    generates a tetrahedral mesh and writes it as an OpenFOAM polyMesh.
    Based on CGAL example mesh_polyhedral_domain_with_features.cpp.

Usage
    \b cgalMeshToPolyMesh inputFile.off

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IOobject.H"
#include "OFstream.H"
#include "meshTools.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "topoSet.H"
#include "processorMeshes.H"

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_criteria_3.h>

// Mesh tools
#include "polyTopoChange.H"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Helper functions

namespace CGALParams = CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

// Convert OpenFOAM triSurface to CGAL Polyhedron directly
void triSurfaceToPolyhedron(const triSurface& surf, Polyhedron& poly)
{
    typedef K::Point_3 Point;
    
    // Get reference to the halfedge data structure
    typedef typename Polyhedron::HalfedgeDS HDS;
    HDS& hds = const_cast<HDS&>(poly.hds());
    
    // Use CGAL's incremental builder to construct the polyhedron directly
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    
    B.begin_surface(surf.nPoints(), surf.size(), 0);
    
    // Add vertices - map OpenFOAM point indices to CGAL vertex handles
    typedef std::map<label, typename HDS::Vertex_handle> VertexMap;
    VertexMap vertexMap;
    
    forAll(surf.points(), i)
    {
        const point& p = surf.points()[i];
        typename HDS::Vertex_handle v = B.add_vertex(Point(p.x(), p.y(), p.z()));
        vertexMap[i] = v;
    }
    
    // Add facets (triangles)
    forAll(surf, facei)
    {
        const triSurface::face_type& f = surf[facei];
        
        if (f.size() != 3)
        {
            FatalErrorInFunction
                << "Expected triangle, got " << f.size() << " vertices"
                << exit(FatalError);
        }
        
        B.begin_facet();
        forAll(f, vi)
        {
            B.add_vertex_to_facet(f[vi]);
        }
        B.end_facet();
    }
    
    B.end_surface();
}

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
    Tr,
    Mesh_domain::Corner_index,
    Mesh_domain::Curve_index> C3t3;

typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

/*
// Convert a CGAL tetrahedral mesh (C3t3) to OpenFOAM polyMesh and write it
void cgalToPolyMesh
(
    const word& regionName,
    const Time& time,
    const fileName& fname
)
{
    // Read input geometry
    Polyhedron polyhedron;

    {
        std::ifstream input(fname);
        if (!input || !input.good())
        {
            FatalErrorInFunction
                << "Cannot open file: " << fname << exit(FatalError);
        }
        input >> polyhedron;

        if (input.fail() || !CGAL::is_triangle_mesh(polyhedron))
        {
            FatalErrorInFunction
                << "Input geometry is not a valid triangle mesh: " << fname
                << exit(FatalError);
        }
    }

    Info<< "Read geometry from " << fname << nl;
    Info<< "Number of vertices: "
        << std::distance(polyhedron.points_begin(), polyhedron.points_end()) << nl;

    // Create domain and detect features
    Mesh_domain domain(polyhedron);
    domain.detect_features();

    Info<< "Features detected" << nl;

    // Mesh criteria (similar to CGAL example)
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    Mesh_criteria criteria
    (
        CGALParams::edge_size(0.025).
        facet_angle(25).facet_size(0.05).facet_distance(0.005).
        cell_radius_edge_ratio(3).cell_size(0.05)
    );

    // Generate mesh
    Info<< "Generating 3D mesh..." << endl;
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    Info<< "Mesh generation done" << nl;

    const Tr& tr = c3t3.triangulation();

    // Map CGAL vertices to point indices
    typedef std::map<typename Tr::Vertex_handle, label> VertexMap;
    VertexMap vertexToIndex;

    pointField points(tr.number_of_vertices());
    label pointi = 0;

    for (typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end(); ++vit)
    {
        const typename K::Point_3& p = vit->point().point();
        points[pointi] = point(p.x(), p.y(), p.z());
        vertexToIndex[vit] = pointi;
        ++pointi;
    }

    // Map CGAL cell handles to cell indices
    typedef std::map<typename Tr::Cell_handle, label> CellMap;
    CellMap cellToIndex;
    for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
         cit != tr.finite_cells_end(); ++cit)
    {
        cellToIndex[cit] = std::distance(tr.finite_cells_begin(), cit);
    }

    // Build faces from tetrahedral cells using finite facets
    DynamicList<face> faces(tr.number_of_cells() * 4);
    DynamicList<label> owner(tr.number_of_cells() * 4);
    DynamicList<label> neighbour(tr.number_of_cells() * 4);

    // Iterate over all finite facets (each geometric face appears exactly once)
    // Use a cell counter instead of std::distance since CGAL iterators don't support it
    label celli = 0;
    for (typename Tr::Finite_facets_iterator fit = tr.finite_facets_begin();
         fit != tr.finite_facets_end(); ++fit, ++celli)
    {
        typename Tr::Cell_handle c1 = fit->first;
        int i1 = fit->second;

        // Get the 3 vertices of this facet from cell c1
        face newFace(3);
        label idx = 0;
        for (label i = 0; i < 4; ++i)
        {
            if (i != i1)
            {
                newFace[idx++] = vertexToIndex[c1->vertex(i)];
            }
        }

        // Check if there's a mirror cell (internal face) or not (boundary face)
        typename Tr::Cell_handle c2 = c1->neighbor(i1);

        bool isInternal = false;
        if (c2 != NULL && !tr.is_infinite(c2))
        {
            isInternal = true;

            // Find the opposite vertex index in c2
            int i2 = -1;
            for (label i = 0; i < 4; ++i)
            {
                if (c2->vertex(i) == c1->vertex(i1))
                {
                    i2 = i;
                    break;
                }
            }

            // Get vertices of the face from c2's perspective
            face mirrorFace(3);
            idx = 0;
            for (label i = 0; i < 4; ++i)
            {
                if (i != i2)
                {
                    mirrorFace[idx++] = vertexToIndex[c2->vertex(i)];
                }
            }

            // Check orientations: in a proper triangulation, the two faces
            // should have opposite orientation for consistent mesh
            bool sameOrientation = (newFace == mirrorFace);

            if (!sameOrientation)
            {
                // Orientations are opposite - this is the expected case
                // Store face with owner=c1, neighbour=c2
                faces.append(newFace);
                owner.append(cellToIndex[c1]);
                neighbour.append(cellToIndex[c2]);
            }
            else
            {
                // Same orientation - swap to ensure consistent convention
                faces.append(mirrorFace);
                owner.append(cellToIndex[c2]);
                neighbour.append(cellToIndex[c1]);
            }
        }

        if (!isInternal)
        {
            // Boundary face (c2 is null, infinite, or not finite)
            faces.append(newFace);
            owner.append(cellToIndex[c1]);
            neighbour.append(-1);
        }
    }

    Info<< "Raw mesh: nFaces=" << faces.size()
        << " nInternal=" << (faces.size() - std::count(neighbour.begin(), neighbour.end(), -1)) << nl;

    // Create polyMesh
    IOobject meshIOObj
    (
        regionName,
        time.constant(),
        time
    );

    polyMesh mesh
    (
        meshIOObj,
        std::move(points),
        std::move(faces.shrink()),
        std::move(owner.shrink()),
        std::move(neighbour.shrink()),
        true  // syncPar
    );

    // Add default boundary patches for any boundary faces
    if (mesh.nFaces() > mesh.nInternalFaces())
    {
        Info<< "Adding default boundary patch" << nl;
        label nBoundaryFaces = mesh.nFaces() - mesh.nInternalFaces();

        polyPatchList patches(1);
        patches.set
        (
            0,
            new polyPatch
            (
                "default",      // name
                nBoundaryFaces, // size
                mesh.nInternalFaces(),  // start
                0,              // index
                mesh.boundaryMesh(),
                word::null,     // physicalType
                wordList(0)     // inGroups
            )
        );
        mesh.addPatches(patches);
    }

    Info<< "Created polyMesh with:" << nl
        << "  nPoints: " << mesh.nPoints() << nl
        << "  nFaces: " << mesh.nFaces() << nl
        << "  nCells: " << mesh.nCells() << nl;

    // Write the mesh
    Info<< nl << "Writing polyMesh to " << time.timePath()/polyMesh::meshSubDir << endl;
    
    // Set write precision
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));
    
    // Remove existing files
    mesh.removeFiles();
    
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed to write polyMesh"
            << exit(FatalError);
    }

    Info<< "Mesh written successfully" << endl;
}
*/

autoPtr<polyMesh> cgalToPolyMesh
(
    const IOobject& meshIOObj,
    const C3t3& c3t3
)
{
    const auto& tr = c3t3.triangulation();

    // Map CGAL vertices to point indices
    typedef std::map<typename C3t3::Triangulation::Vertex_handle, label> VertexMap;
    VertexMap vertexToIndex;

    pointField points(tr.number_of_vertices());
    label pointi = 0;

    for (typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end(); ++vit)
    {
        const typename K::Point_3& p = vit->point().point();
        points[pointi] = point(p.x(), p.y(), p.z());
        vertexToIndex[vit] = pointi;
        ++pointi;
    }
    Info<< "points:" << pointi << endl;

    // Map CGAL cell handles to cell indices
    const cellModel& tet = cellModel::ref(cellModel::TET);
    labelList tetPoints(4);

    typedef typename C3t3::Cells_in_complex_iterator Cell_iterator;
    
    cellShapeList cells(c3t3.number_of_cells());
    Info<< "cells:" << cells.size() << endl;
    for
    (
        Cell_iterator cit = c3t3.cells_in_complex_begin();
         cit != c3t3.cells_in_complex_end(); ++cit)
    {
        const label celli =
            std::distance(c3t3.cells_in_complex_begin(), cit);
        // cellToIndex[cit] = celli;

        typename Tr::Vertex_handle v0 = cit->vertex(0);
        typename Tr::Vertex_handle v1 = cit->vertex(1);
        typename Tr::Vertex_handle v2 = cit->vertex(2);
        typename Tr::Vertex_handle v3 = cit->vertex(3);
    
        tetPoints[0] = vertexToIndex[v0];
        tetPoints[1] = vertexToIndex[v1];
        tetPoints[2] = vertexToIndex[v2];
        tetPoints[3] = vertexToIndex[v3];
        cells[celli].reset(tet, tetPoints);
        // const typename K::Point_3& p0 = v0->point().point();
        // const typename K::Point_3& p1 = v1->point().point();
        // const typename K::Point_3& p2 = v2->point().point();
        // const typename K::Point_3& p3 = v3->point().point();
        // Pout<< "cell:" << celli << " verts:" << cells[celli] << endl;
        // celli++;
    }
    // cells.setSize(celli);

    // Create polyMesh
    return autoPtr<polyMesh>::New
    (
        meshIOObj,
        std::move(points),
        cells,
        faceListList(),
        wordList(),
        wordList(),
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Read geometry using CGAL, extract features, generate tetrahedral mesh"
    );

    argList::addArgument("inputFile", "The input geometry file (OFF format)");

    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"
    #include "createTime.H"

    const fileName importName = args.get<fileName>(1);

    Info<< "Converting CGAL mesh to OpenFOAM polyMesh" << nl;
    Info<< "Input file: " << importName << nl;

    const word oldInstance = runTime.instance();
    DebugVar(oldInstance);


    // // Create constant directory
    // fileName constantDir = runTime.path()/"constant";
    // fileName meshDir = constantDir/polyMesh::meshSubDir;



   // Read input geometry
   // ~~~~~~~~~~~~~~~~~~~

    const triSurfaceMesh surf
    (
        IOobject
        (
            importName,         // name
            runTime.constant(),   // instance
            "triSurface",         // local
            runTime,              // registry
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );
    
    // Convert triSurface to Polyhedron for CGAL meshing
    Polyhedron polyhedron;
    triSurfaceToPolyhedron(surf, polyhedron);


    Info<< "Read geometry from " << importName << nl;
    Info<< "Number of vertices: "
        << std::distance(polyhedron.points_begin(), polyhedron.points_end())
        << nl;


    // Generate mesh
    // ~~~~~~~~~~~~~

    // Create domain and detect features
    Mesh_domain domain(polyhedron);
    domain.detect_features();

    Info<< "Features detected" << nl;

    // Mesh criteria (similar to CGAL example)
    
    Mesh_criteria criteria
    (
        CGALParams::edge_size(0.025).
        facet_angle(25).facet_size(0.05).facet_distance(0.005).
        cell_radius_edge_ratio(3).cell_size(0.05)
    );

    // Generate mesh
    Info<< "Generating 3D mesh..." << endl;
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    Info<< "Mesh generation done" << nl;



    // Convert to polyMesh
    autoPtr<polyMesh> meshPtr
    (
        cgalToPolyMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime.constant(),
                runTime
            ),
            c3t3
        )
    );
    auto& mesh = *meshPtr;


    Info<< "Finding correspondence to surface" << nl;

    pointField bfc
    (
        SubList<point>
        (
            mesh.faceCentres(),
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        )
    );
    List<pointIndexHit> info;
    surf.findNearest
    (
        bfc,
        scalarField(bfc.size(), mesh.bounds().magSqr()),
        info
    );
    labelList surfaceFaceIDs(mesh.nBoundaryFaces(), -1);
    // Per boundary faces the patch. Start from patch0 (=defaultFaces)
    labelList surfaceRegionIDs(mesh.nBoundaryFaces(), 0);
    forAll(info, i)
    {
        if (info[i].hit())
        {
            surfaceFaceIDs[i] = info[i].index();
            const labelledTri& tri =
                static_cast<const triSurface&>(surf)[surfaceFaceIDs[i]];
            surfaceRegionIDs[i] = tri.region();
        }
    }



    // Add zero-sized patches
    // ~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Adding zero-sized patches" << nl;

    const auto& surfPatches = surf.patches();
    polyPatchList patches(surfPatches.size());
    label startFace = mesh.nInternalFaces();
    label nFaces = mesh.nBoundaryFaces();
    forAll(surfPatches, patchi)
    {
        const auto& surfPatch = surfPatches[patchi];

        Info<< "    adding patch " << surfPatch.name() << nl;
        patches.set
        (
            patchi,
            new polyPatch
            (
                surfPatch.name(),
                nFaces,
                startFace,
                patchi,
                mesh.boundaryMesh(),
                surfPatch.geometricType()
            )
        );
        startFace += nFaces;
        nFaces = 0;
    }
    mesh.removeBoundary();
    mesh.addPatches(patches);



    // Move boundary faces to correct patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Moving boundary faces to correct patches" << nl;

    polyTopoChange meshMod(mesh);

    DynamicList<label> zones;
    DynamicList<bool> flips;
    forAll(surfaceRegionIDs, i)
    {
        const label facei = i + mesh.nInternalFaces();
        const label patchi = surfaceRegionIDs[i];
        meshMod.faceZones(facei, zones, flips);

        // Modify the face with new patch assignment
        meshMod.modifyFace
        (
            mesh.faces()[facei],
            facei,
            mesh.faceOwner()[facei],
            -1,
            false,  // flipFaceFlux
            patchi,
            zones,
            flips
        );
    }


    // Create mesh, return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map());

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    // mesh.setInstance(oldInstance);
    mesh.setInstance(runTime.constant());

    Info<< "Writing mesh to " << mesh.pointsInstance() << endl;

    mesh.write();
    topoSet::removeFiles(mesh);
    processorMeshes::removeFiles(mesh);

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
