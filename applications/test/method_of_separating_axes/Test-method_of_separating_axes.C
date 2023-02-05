/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Implementation of
    "Intersection of Convex Objects: The Method of Separating Axes",
    David Eberly, Geometric Tools, LLC

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "polyTopoChange.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "columnFvMesh.H"
#include "OBJstream.H"
#include "regionSplit.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int WhichSide
(
    const pointField& points,
    const labelUList& S,
    const vector& D,
    const point& P
)
{
    // S vertices are projected to the form P+t*D. Return value is +1 if all
    // t > 0, -1 if all t < 0, 0 otherwise, in which case the line splits
    // the polygon.

    int positive = 0;
    int negative = 0;

    forAll(S, i)
    {
        const scalar t = (D & (points[S[i]]-P));
        if (t > 0)
        {
            positive++;
        }
        else if (t < 0)
        {
            negative++;
        }
        if (positive && negative)
        {
            return 0;
        }
    }
    return (positive ? +1 : -1);
}


bool TestIntersection(const polyMesh& mesh, const label C0, const label C1)
{
    const auto& cells = mesh.cells();
    //const auto& faces = mesh.faces();
    const auto& points = mesh.points();
    const auto& faceAreas = mesh.faceAreas();
    const auto& faceOwner = mesh.faceOwner();
    const auto& edges = mesh.edges();

    // Test faces of C0 for separation. Because of the counterclockwise
    // ordering, the projection interval for C0 is [m,0] where m <= 0.
    // Only try to determine if C1 is on the 'positive' side of the line.
    for (const label facei : cells[C0])
    {
        //const face& f = faces[facei];
        const vector D
        (
            C0 == faceOwner[facei]
          ? faceAreas[facei]
          : -faceAreas[facei]
        );
        const point& fc = mesh.faceCentres()[facei];

        if (WhichSide(points, mesh.cellPoints(C1), D, fc) > 0)
        {
            // C1 is entirely on 'positive' side of line
            Pout<< "cell:" << C1 << " bb:"
                << boundBox(points, mesh.cellPoints(C1))
                << " is completely positive of plane through fc:" << fc
                << " direction:" << D/mag(D) << endl;

            return false;
        }
    }

    // Test faces of C1 for separation. Because of the counterclockwise
    // ordering, the projection interval for C1 is [m,0] where m <= 0.
    // Only try to determine if C0 is on the 'positive' side of the line.
    for (const label facei : cells[C1])
    {
        //const face& f = faces[facei];
        const vector D
        (
            C1 == faceOwner[facei]
          ? faceAreas[facei]
          : -faceAreas[facei]
        );
        const point& fc = mesh.faceCentres()[facei];
        if (WhichSide(points, mesh.cellPoints(C0), D, fc) > 0)
        {
            // C0 is entirely on 'positive' side of line
            Pout<< "cell:" << C0 << " bb:"
                << boundBox(points, mesh.cellPoints(C0))
                << " is completely positive of plane through fc:" << fc
                << " direction:" << D/mag(D) << endl;
            return false;
        }
    }

    // Test cross product of pairs of edges, one from each polyhedron
    for (const label edge0 : mesh.cellEdges(C0))
    {
        const edge& e0 = edges[edge0];
        const point& e0Pt = points[e0.start()];

        for (const label edge1 : mesh.cellEdges(C1))
        {
            const edge& e1 = edges[edge1];
            const vector D(e0.vec(points) ^ e1.vec(points));
            int side0 = WhichSide(points, mesh.cellPoints(C0), D, e0Pt);
            if (side0 == 0)
            {
                continue;
            }

            int side1 = WhichSide(points, mesh.cellPoints(C1), D, e0Pt);
            if (side1 == 0)
            {
                continue;
            }

            if (side0*side1 < 0)
            {
                 // C0 and C1 are on 'opposite' sides of line e0Pt+t*D
                return false;
            }
        }
    }
    return true;
}


int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    OBJstream os0(mesh.time().path()/"intersects0.obj");
    Pout<< "Dumping intersecting cells to " << os0.name() << endl;
    OBJstream os1(mesh.time().path()/"intersects1.obj");
    Pout<< "Dumping intersecting cells to " << os1.name() << endl;
    OBJstream osLines(mesh.time().path()/"intersectionLines.obj");
    Pout<< "Dumping intersecting cells to " << osLines.name() << endl;

    const auto& ccs = mesh.cellCentres();
    const auto& regions = regionSplit::New(mesh);

    for (label r0 = 0; r0 < regions.nRegions()-1; r0++)
    {
        forAll(ccs, C0)
        {
            if (regions[C0] == r0)
            {
                Pout<< "Region0 : C0:" << C0 << " cc:" << ccs[C0] << endl;
                for (label r1 = r0+1; r1 < regions.nRegions(); r1++)
                {
                    forAll(ccs, C1)
                    {
                        if (regions[C1] == r1)
                        {
                            Pout<< "Region1 : C1:" << C1 << " cc:" << ccs[C1]
                                << endl;
                            bool intersects = TestIntersection(mesh, C0, C1);
                            Pout<< "    intersect:" << intersects << endl;

                            if (intersects)
                            {
                                const auto& faces = mesh.faces();
                                const auto& points = mesh.points();
                                const cell& C0Faces = mesh.cells()[C0];
                                os0.write(UIndirectList<face>(faces, C0Faces)(), points, false);
                                const cell& C1Faces = mesh.cells()[C1];
                                os1.write(UIndirectList<face>(faces, C1Faces)(), points, false);
                                osLines.write(linePointRef(ccs[C0], ccs[C1]));
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}

// ************************************************************************* //
