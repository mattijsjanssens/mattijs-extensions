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

Application
    Test-lduSolver

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "lduMatrix.H"
#include "IFstream.H"
#include "DynamicField.H"
#include "lduPrimitiveMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    //#include "createMesh.H"

    IFstream is("matrixmarket.txt");
    if (!is.good())
    {
        FatalIOErrorInFunction(is) << "Problem opening file"
            << exit(FatalIOError);
    }
    // %%MatrixMarket
    word hdr(is);
    DebugVar(hdr);

    // matrix
    word type(is);
    DebugVar(type);

    // coordinate
    word coord(is);
    DebugVar(coord);

    // real
    word numType(is);
    DebugVar(numType);

    // general
    word matrixType(is);
    DebugVar(matrixType);

    const label nRows = readLabel(is);
    const label nCols = readLabel(is);
    if (nRows != nCols)
    {
        FatalIOErrorInFunction(is) << "Non-square matrix : rows:" << nRows
            << " columns:" << nCols << exit(FatalIOError);
    }
    const label nEntries = readLabel(is);
    DebugVar(nEntries);

    // Read
    DynamicList<label> row(nEntries);
    DynamicList<label> column(nEntries);
    DynamicField<scalar> coeff(nEntries);

    for (label conni = 0; conni < nEntries; conni++)
    {
        row.append(readLabel(is)-1);
        column.append(readLabel(is)-1);
        coeff.append(readScalar(is));
    }

    // Count nNbrs per row
    labelList nLower(nRows, 0);
    labelList nUpper(nRows, 0);

    // And extract the diagonal since we do the checks anyway
    scalarField diag(nRows, 0.0);

    forAll(row, conni)
    {
        const label celli = row[conni];
        if (column[conni] < celli)
        {
            nLower[celli]++;
        }
        else if (column[conni] > celli)
        {
            nUpper[celli]++;
        }
        else
        {
            diag[celli] = coeff[conni];
        }
    }

DebugVar(nLower);
DebugVar(nUpper);
DebugVar(diag);



    // Extract out the connections to higher numbered cells
    //const label nFaces = (nEntries - nRows)/2;
    const label nFaces = max
    (
        sum(nUpper),
        sum(nLower)
    );
DebugVar(nFaces);


    DynamicList<label> lowerAddr(nFaces);
    DynamicField<scalar> lower(nFaces);
    DynamicList<label> upperAddr(nFaces);
    DynamicList<scalar> upper(nFaces);


//     // Collect all higher numbered cells
//     labelListList upperCells(nRows);
//     forAll(upperCells, celli)
//     {
//         upperCells[celli].setSize(nUpper[celli]);
//     }
//     nUpper = 0;
//     forAll(row, conni)
//     {
//         label celli = row[conni];
//         label nbri = column[conni];
//         if (nbri > celli)
//         {
//             upperCells[celli][nUpper[celli]++] = nbri;
//         }
//     }
// 
//     // Allocate faces to higher numbered cells
//     forAll(upperCells, celli)
//     {
//         labelList& u = upperCells[celli];
//         sort(u);
//         forAll(u, i)
//         {
//             // Allocate a face
//             label facei = lowerAddr.size();
//             lowerAddr.append(celli);
//             upperAddr.append(u[i]);
//             upper.append(coeff[conni]);
//             lower.append(-1);
// 





    labelListList lowerCells(nRows);
    forAll(lowerCells, celli)
    {
        lowerCells[celli].setSize(nLower[celli]);
    }
    nLower = 0;

    // Allocate faces to higher numbered cells (assumes are output in increasing
    // column number
    forAll(row, conni)
    {
        label celli = row[conni];
        label nbri = column[conni];
        if (nbri > celli)
        {
            // Allocate a face
            label facei = lowerAddr.size();
            lowerAddr.append(celli);
            upperAddr.append(nbri);
            upper.append(coeff[conni]);
            lower.append(-1);

            // Store the face
            lowerCells[nbri][nLower[nbri]++] = facei;
        }
    }

    // Use the addressing to set the coefficients to lower numbered cells
    forAll(row, conni)
    {
        label celli = row[conni];
        label nbri = column[conni];
        if (nbri < celli)
        {
            // Find the face
            const labelList& faces = lowerCells[celli];

            label facei = -1;
            forAll(faces, i)
            {
                if (lowerAddr[faces[i]] == nbri)
                {
                    facei = faces[i];
                    break;
                }
            }

            if (facei == -1)
            {
                FatalErrorInFunction << "Asymetric addressing not supported."
                    << " Cannot find face between " << celli
                    << " and " << nbri << exit(FatalError);
            }

            lower[facei] = coeff[conni];
        }
    }

    DebugVar(lowerAddr);
    DebugVar(lower);
    DebugVar(upperAddr);
    DebugVar(upper);


    PtrList<const lduInterface> interfaces(0);
    lduPrimitiveMesh lMesh
    (
        nRows,
        lowerAddr,
        upperAddr,
        interfaces,
        lduSchedule(0),
        Pstream::worldComm
    );

    lduMatrix lMat(lMesh);
    lMat.diag() = diag;
    lMat.lower() = lower;
    lMat.upper() = upper;

    // No boundary handling
    const lduInterfaceFieldPtrsList scalarInterfaces(0);
    const FieldField<Field, scalar> interfaceBouCoeffs(0);
    const FieldField<Field, scalar> interfaceIntCoeffs(0);

    // Source
    const scalarField totalSource(nRows, 0.0);

    // Solution
    scalarField psi(nRows, 123.0);

    dictionary solverControls;

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        "p",
        lMat,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        scalarInterfaces,
        solverControls
    )->solve(psi, totalSource);

    solverPerf.print(Info.masterStream(lMesh.comm()));

    return 0;
}


// ************************************************************************* //
