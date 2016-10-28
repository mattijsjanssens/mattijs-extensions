/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>

using std::ofstream;
using std::ios;

#include "Time.H"
#include "fluentFvMesh.H"
#include "primitiveMesh.H"
#include "wallFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "symmetryFvPatch.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluentFvMesh::fluentFvMesh(const IOobject& io)
:
    fvMesh(io)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluentFvMesh::writeCellShapes
(
    ostream& fluentMeshFile,
    const cellShapeList& cells
) const
{
    static bool hasWarned = false;

    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& tet = *(cellModeller::lookup("tet"));

    forAll(cells, celli)
    {
        if (cells[celli].model() == tet)
        {
            fluentMeshFile << " " << 2;
        }
        else if (cells[celli].model() == hex)
        {
            fluentMeshFile << " " << 4;
        }
        else if (cells[celli].model() == pyr)
        {
            fluentMeshFile << " " << 5;
        }
        else if (cells[celli].model() == prism)
        {
            fluentMeshFile << " " << 6;
        }
        else
        {
            if (!hasWarned)
            {
                hasWarned = true;

                WarningInFunction
                    << "foamMeshToFluent: cell shape for cell "
                    << celli << " only supported by Fluent polyhedral meshes."
                    << nl
                    << "    Suppressing any further messages for polyhedral"
                    << " cells." << endl;
            }
            fluentMeshFile << " " << 7;
        }
    }
}


void Foam::fluentFvMesh::writeFluentMesh() const
{
    // make a directory called proInterface in the case
    mkDir(time().rootPath()/time().caseName()/"fluentInterface");

    // open a file for the mesh
    ofstream fluentMeshFile
    (
        (
            time().rootPath()/
            time().caseName()/
            "fluentInterface"/
            time().caseName() + ".msh"
        ).c_str()
    );

    Info<< "Writing Header" << endl;

    fluentMeshFile
        << "(0 \"FOAM to Fluent Mesh File\")" << std::endl << std::endl
        << "(0 \"Dimension:\")" << std::endl
        << "(2 3)" << std::endl << std::endl
        << "(0 \"Grid dimensions:\")" << std::endl;

    // Writing number of points
    fluentMeshFile
            << "(10 (0 1 ";

    // Writing hex
    fluentMeshFile.setf(ios::hex, ios::basefield);

    fluentMeshFile
        << nPoints() << " 0 3))" << std::endl;

    // Writing number of cells
    fluentMeshFile
        << "(12 (0 1 "
        << nCells() << " 0 0))" << std::endl;

    // Writing number of faces
    label nFcs = nFaces();

    fluentMeshFile
            << "(13 (0 1 ";

    // Still writing hex
    fluentMeshFile
        << nFcs << " 0 0))" << std::endl << std::endl;

    // Return to dec
    fluentMeshFile.setf(ios::dec, ios::basefield);

    // Writing points
    fluentMeshFile
            << "(10 (1 1 ";

    fluentMeshFile.setf(ios::hex, ios::basefield);
    fluentMeshFile
        << nPoints() << " 1 3)"
        << std::endl << "(" << std::endl;

    fluentMeshFile.precision(10);
    fluentMeshFile.setf(ios::scientific);

    const pointField& p = points();

    forAll(p, pointi)
    {
        fluentMeshFile
            << "    "
            << p[pointi].x() << " "
            << p[pointi].y()
            << " " << p[pointi].z() << std::endl;
    }

    fluentMeshFile
        << "))" << std::endl << std::endl;

    const labelUList& own = owner();
    const labelUList& nei = neighbour();

    const faceList& fcs = faces();

    // Writing (mixed) internal faces
    fluentMeshFile
        << "(13 (9 1 "
        << own.size() << " 2 0)" << std::endl << "(" << std::endl;

    forAll(own, facei)
    {
        const labelList& l = fcs[facei];

        fluentMeshFile << "    ";

        fluentMeshFile << l.size() << " ";

        forAll(l, lI)
        {
            fluentMeshFile << l[lI] + 1 << " ";
        }

        fluentMeshFile << nei[facei] + 1 << " ";
        fluentMeshFile << own[facei] + 1 << std::endl;
    }

    fluentMeshFile << "))" << std::endl;

    label nWrittenFaces = own.size();

    // Writing boundary faces
    forAll(boundary(), patchi)
    {
        const faceUList& patchFaces = boundaryMesh()[patchi];

        const labelList& patchFaceCells =
            boundaryMesh()[patchi].faceCells();

        // The face group will be offset by 10 from the patch label

        // Write header
        fluentMeshFile
            << "(13 (" << patchi + 10 << " " << nWrittenFaces + 1
            << " " << nWrittenFaces + patchFaces.size() << " ";

        nWrittenFaces += patchFaces.size();

        // Write patch type
        if (isA<wallFvPatch>(boundary()[patchi]))
        {
            fluentMeshFile << 3;
        }
        else if
        (
            isA<symmetryPlaneFvPatch>(boundary()[patchi])
         || isA<symmetryFvPatch>(boundary()[patchi])
        )
        {
            fluentMeshFile << 7;
        }
        else
        {
            fluentMeshFile << 4;
        }

        fluentMeshFile
            <<" 0)" << std::endl << "(" << std::endl;

        forAll(patchFaces, facei)
        {
            const labelList& l = patchFaces[facei];

            fluentMeshFile << "    ";

            fluentMeshFile << l.size() << " ";

            // Note: In Fluent, all boundary faces point inwards, which is
            // opposite from the OpenFOAM convention.
            // Turn them around on printout
            forAllReverse (l, lI)
            {
                fluentMeshFile << l[lI] + 1 << " ";
            }

            fluentMeshFile << patchFaceCells[facei] + 1 << " 0" << std::endl;
        }

        fluentMeshFile << "))" << std::endl;
    }


    // Determine if cellZones and ordered.
    const cellZoneMesh& czm = cellZones();

    bool ordered;
    labelList nZoneCells(czm.size(), 0);
    {
        labelList zoneID(nCells(), -1);

        forAll(czm, zonei)
        {
            UIndirectList<label>(zoneID, czm[zonei]) = zonei;
            nZoneCells[zonei] = czm[zonei].size();
        }

        // Sort in increasing order
        labelList newToOld;
        sortedOrder(zoneID, newToOld);

        // Check for consistent order
        ordered = true;
        label prevZonei = zoneID[0];
        for (label celli = 1; celli < zoneID.size(); celli++)
        {
            if (zoneID[celli] < prevZonei)
            {
                ordered = false;
                break;
            }
            prevZonei = zoneID[celli];
        }

        if (ordered)
        {
            Info<< "Detected mesh with ordered cellZones."
                << "Writing Fluent mesh output with cellZones." << endl;
        }
        else if (czm.size())
        {
            Info<< "Detected mesh with unordered cellZones."
                << " Ignoring zone info. Run renumberMesh "
                << "with cellZoneOrder method to enable Fluent mesh"
                << " output with cellZones." << endl;

            nZoneCells.clear();
        }
    }

DebugVar(ordered);
DebugVar(nZoneCells);



    // Writing cells
    // ~~~~~~~~~~~~~

    label nUnzonedCells = nCells()-sum(nZoneCells);

    // Write unzoned cells first
    {
        fluentMeshFile
            << "(12 (1 1 "
            << nUnzonedCells << " 1 0)(" << std::endl;

        SubList<cellShape> set(cellShapes(), nUnzonedCells);
        writeCellShapes(fluentMeshFile, set);

        fluentMeshFile << ")())" << std::endl;
        fluentMeshFile << "(39 (1 fluid fluid-1)())" << std::endl;
    }

    label startCelli = nUnzonedCells;
    forAll(nZoneCells, zonei)
    {
        fluentMeshFile
            << "(12 (" << zonei+1 <<  " 1 "
            << nZoneCells[zonei] << " 1 0)(" << std::endl;

        SubList<cellShape> set(cellShapes(), nZoneCells[zonei], startCelli);
        writeCellShapes(fluentMeshFile, set);

        fluentMeshFile << ")())" << std::endl;
        fluentMeshFile << "(39 (" << zonei+1 << " fluid " << czm[zonei].name()
            << ")())" << std::endl;

        startCelli += nZoneCells[zonei];
    }


    // Return to dec
    fluentMeshFile.setf(ios::dec, ios::basefield);

    // Internal faces
    fluentMeshFile << "(39 (9 interior interior-1)())" << std::endl;

    // Writing boundary patch types
    forAll(boundary(), patchi)
    {
        fluentMeshFile
            << "(39 (" << patchi + 10 << " ";

        // Write patch type
        if (isA<wallFvPatch>(boundary()[patchi]))
        {
            fluentMeshFile << "wall ";
        }
        else if
        (
            isA<symmetryPlaneFvPatch>(boundary()[patchi])
         || isA<symmetryFvPatch>(boundary()[patchi])
        )
        {
            fluentMeshFile << "symmetry ";
        }
        else
        {
            fluentMeshFile << "pressure-outlet ";
        }

        fluentMeshFile
            << boundary()[patchi].name() << ")())" << std::endl;
    }
}


// ************************************************************************* //
