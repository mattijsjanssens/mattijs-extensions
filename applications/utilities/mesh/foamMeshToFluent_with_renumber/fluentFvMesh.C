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
#include "cellZoneRenumber.H"

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
    // Layout of Fluent zones:
    //  1               : unzoned cells
    //  2..zonei        : cellZones
    //  2+nZones        : interior faces
    //  3+nZones+patchi : patch faces

    const cellZoneMesh& czm = cellZones();


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
        << "(13 (" << 2+czm.size() << " 1 "
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

        // The face group will be offset by cellzones from the patch label

        // Write header
        fluentMeshFile
            << "(13 (" << 3+czm.size()+patchi << " " << nWrittenFaces + 1
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
                << "with the " << cellZoneRenumber::typeName
                << " renumber method to enable Fluent mesh"
                << " output with cellZones." << endl;

            nZoneCells.clear();
        }
    }



    // Writing cells
    // ~~~~~~~~~~~~~

    label nUnzonedCells = nCells()-sum(nZoneCells);


    label nWrittenCells = 0;

    // Write unzoned cells first
    {
        fluentMeshFile
            << "(12 (1 1 "
            << nUnzonedCells << " 1 0)(" << std::endl;

        SubList<cellShape> set(cellShapes(), nUnzonedCells);
        writeCellShapes(fluentMeshFile, set);

        nWrittenCells += nUnzonedCells;

        fluentMeshFile << ")())" << std::endl;
        fluentMeshFile << "(39 (1 fluid fluid-1)())" << std::endl;
    }

    forAll(nZoneCells, zonei)
    {
        fluentMeshFile
            << "(12 (" << 2+zonei << " " << nWrittenCells+1 << " "
            << nWrittenCells+nZoneCells[zonei] << " 1 0)(" << std::endl;

        SubList<cellShape> set(cellShapes(), nZoneCells[zonei], nWrittenCells);
        writeCellShapes(fluentMeshFile, set);

        fluentMeshFile << ")())" << std::endl;
        fluentMeshFile << "(39 (" << 2+zonei << " fluid " << czm[zonei].name()
            << ")())" << std::endl;

        nWrittenCells += nZoneCells[zonei];
    }


    // Return to dec
    fluentMeshFile.setf(ios::dec, ios::basefield);

    // Internal faces
    fluentMeshFile << "(39 (" << 2+czm.size() << " interior interior-1)())"
        << std::endl;

    // Writing boundary patch types
    forAll(boundary(), patchi)
    {
        fluentMeshFile
            << "(39 (" << 3+czm.size()+patchi << " ";

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
