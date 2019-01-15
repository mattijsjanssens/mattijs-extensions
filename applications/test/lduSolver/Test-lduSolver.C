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
#include "CompactListList.H"
#include "scalarIOField.H"
#include "IOdictionary.H"
#include "lduPrimitiveInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void readMatrixMarket
(
    const fileName& fName,

    label& nRows,
    labelList& lowerAddr,
    labelList& upperAddr,
    PtrList<const lduInterface>& interfaces,

    scalarField& diag,
    scalarField& lower,
    scalarField& upper
)
{
    IFstream is(fName);
    if (!is.good())
    {
        FatalIOErrorInFunction(is) << "Problem opening " << fName
            << exit(FatalIOError);
    }
    // %%MatrixMarket
    word hdr(is);
    DebugVar(hdr);
    if (hdr != "%%MatrixMarket")
    {
        FatalIOErrorInFunction(is) << "Incorrect header " << hdr
            << exit(FatalIOError);
    }

    // matrix
    word object(is);
    DebugVar(object);
    if (object != "matrix")
    {
        FatalIOErrorInFunction(is) << "Unsupported object type " << object
            << exit(FatalIOError);
    }

    // coordinate
    word format(is);
    DebugVar(format);
    if (format != "coordinate")
    {
        FatalIOErrorInFunction(is) << "Unsupported format " << format
            << exit(FatalIOError);
    }

    // real
    word fieldType(is);
    DebugVar(fieldType);
    if (fieldType != "real")
    {
        FatalIOErrorInFunction(is) << "Unsupported field " << fieldType
            << exit(FatalIOError);
    }

    // general
    word symmetryType(is);
    DebugVar(symmetryType);
    if (symmetryType != "general")
    {
        FatalIOErrorInFunction(is) << "Format " << symmetryType
            << " not supported. Only supported format is 'general'"
            << exit(FatalIOError);
    }


    nRows = readLabel(is);
    const label nCols = readLabel(is);
    if (nRows != nCols)
    {
        FatalIOErrorInFunction(is) << "Non-square matrix : rows:" << nRows
            << " columns:" << nCols << exit(FatalIOError);
    }
    const label nEntries = readLabel(is);
    DebugVar(nEntries);

    // Read into structure
    typedef Tuple2<labelPair, scalar> CellsCoeffType;
    DynamicList<CellsCoeffType> entries(nEntries);
    for (label conni = 0; conni < nEntries; conni++)
    {
        const label row = readLabel(is)-1;
        const label col = readLabel(is)-1;
        const scalar coeff = readScalar(is);
        entries.append(CellsCoeffType(labelPair(row, col), coeff));
    }

    // Do ordering into increasing nbr so when we hand out faces they
    // are in upper-triangular order
    Foam::sort
    (
        entries,
        [](const CellsCoeffType& i, const CellsCoeffType& j)
        {
            const labelPair& ei = i.first();
            const labelPair& ej = j.first();

            if (ei.first() < ej.first())
            {
                return true;
            }
            else if (ei.first() == ej.first())
            {
                if (ei.second() == ej.second())
                {
                    FatalErrorInFunction << "Problem : multiple coefficients"
                        << " between cell " << ei.first() << " and "
                        <<  ei.second() << exit(FatalError);
                }
                return ei.second() < ej.second();
            }
            else
            {
                return false;
            }
        }
    );

    // Count nNbrs per row
    labelList nLower(nRows, 0);
    labelList nUpper(nRows, 0);

    // And extract the diagonal since we do the checks anyway
    diag.setSize(nRows);
    diag = 0.0;

    forAll(entries, conni)
    {
        const CellsCoeffType& e = entries[conni];
        const label celli = e.first().first();
        if (e.first().second() < celli)
        {
            nLower[celli]++;
        }
        else if (e.first().second() > celli)
        {
            nUpper[celli]++;
        }
        else
        {
            diag[celli] = e.second();
        }
    }

    // Extract out the connections to higher numbered cells
    //const label nFaces = (nEntries - nRows)/2;
    const label nFaces = max
    (
        sum(nUpper),
        sum(nLower)
    );

    label facei = 0;
    lowerAddr.setSize(nFaces);
    lower.setSize(nFaces);
    upperAddr.setSize(nFaces);
    upper.setSize(nFaces);

    // Allocate faces to higher numbered cells and store correspondence

    CompactListList<label> lowerCells(nLower);
    const labelList& offsets = lowerCells.offsets();
    labelList& m = lowerCells.m();

    nLower = 0;
    forAll(entries, conni)
    {
        const CellsCoeffType& e = entries[conni];

        label celli = e.first().first();
        label nbri = e.first().second();
        if (nbri > celli)
        {
            // Allocate a face
            lowerAddr[facei] = celli;
            upperAddr[facei] = nbri;
            upper[facei] = e.second();
            lower[facei] = -1;
            // Store for later lookup
            m[offsets[nbri]+nLower[nbri]++] = facei;

            facei++;
        }
    }


    // Use the addressing to set the coefficients to lower numbered cells
    forAll(entries, conni)
    {
        const CellsCoeffType& e = entries[conni];

        label celli = e.first().first();
        label nbri = e.first().second();
        if (nbri < celli)
        {
            // Find the face
            const labelUList& faces = lowerCells[celli];

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

            lower[facei] = e.second();
        }
    }

//DebugVar(lowerAddr);
//DebugVar(lower);
//DebugVar(upperAddr);
//DebugVar(upper);
}

int main(int argc, char *argv[])
{
    argList::validArgs.append("matrixmarket");
    argList::validArgs.append("source");
    #include "setRootCase.H"
    #include "createTime.H"

    fileName fName(args.args()[1]);
    fileName sourceName(args.args()[2]);

    Info<< "Reading matrix from " << fName << nl
        << "        source from " << sourceName << nl << endl;

    label nCells;
    labelList lowerAddr;
    labelList upperAddr;
    PtrList<const lduInterface> interfaces;

    scalarField diag;
    scalarField lower;
    scalarField upper;

    // Read Matrix market
    {
        readMatrixMarket
        (
            fName,

            // Addressing
            nCells,
            lowerAddr,
            upperAddr,
            interfaces,

            // Coefficients
            diag,
            lower,
            upper
        );
    }

    
    // Source
    scalarField totalSource;
    {
        const IOdictionary sourceDict
        (
            IOobject
            (
                sourceName,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        totalSource = scalarField("value", sourceDict, nCells);
    }


    const lduPrimitiveMesh lMesh
    (
        nCells,
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

    // Solution
    scalarField psi(nCells, 123.0);

    dictionary solverControls("Test-lduSolver");
    solverControls.add("solver", "PBiCGStab");
    solverControls.add("preconditioner", "DILU");
    solverControls.add("tolerance", 1e-06);
    solverControls.add("relTol", 0);

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

    DebugVar(psi);



    if (Pstream::parRun())
    {
        labelList decomposition;
        if (Pstream::master())
        {
            // Decompose in identical parts
            const label nLocal = nCells / Pstream::nProcs();
            labelList offsets(Pstream::nProcs()+1);
            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                offsets[proci+1] = offsets[proci]+nLocal;
            }
            // Any remainder into last processor
            offsets.last() = nCells;

            decomposition.setSize(nCells);
            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                label size = offsets[proci+1]-offsets[proci];
                SubList<label>(decomposition, size, offsets[proci]) = proci;
            }
        }

        // Assign local cell numbers



        //const globalIndex globalNumbering(offsets);

        // From my processor cell to global cell
        List<DynamicList<label>> procCellMap(Pstream::nProcs());
        // Local cell number
        labelList localCell(decomposition.size());
        forAll(decomposition, celli)
        {
            localCell[celli] = procCellMap[decomposition[celli]].size();
            procCellMap[decomposition[celli]].append(celli);
        }

        // From my processor to global face
        List<DynamicList<label>> procFaceMap(Pstream::nProcs());

        // From my processor and nbr processor to global faces
        List<List<DynamicList<label>>> procNbrFaces(Pstream::nProcs());
        List<List<DynamicList<label>>> procFaceCells(Pstream::nProcs());
        forAll(procNbrFaces, proci)
        {
            procNbrFaces[proci].setSize(Pstream::nProcs());
            procFaceCells[proci].setSize(Pstream::nProcs());
        }

        forAll(lMesh.lowerAddr(), facei)
        {
            label ownCelli = lMesh.lowerAddr()[facei];
            label ownProc = decomposition[ownCelli];
            label ownProcCelli = localCell[ownCelli];
            label nbrCelli = lMesh.upperAddr()[facei];
            label nbrProc = decomposition[nbrCelli];
            label nbrProcCelli = localCell[nbrCelli];

            if (ownProc == nbrProc)
            {
                // Internal face
                procFaceMap[ownProc].append(facei);
            }
            else
            {
                // Inter-processor map
                procNbrFaces[ownProc][nbrProc].append(facei);
                procFaceCells[ownProc][nbrProc].append(ownProcCelli);
                procNbrFaces[nbrProc][ownProc].append(facei);
                procFaceCells[nbrProc][ownProc].append(nbrProcCelli);
            }
        }

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (Pstream::master())
        {
            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                UOPstream os(proci, pBufs);

                os  << nCells
                    <<  UIndirectList<label>
                        (
                            lMesh.lowerAddr(),
                            procFaceMap[proci]
                        )
                    <<  UIndirectList<label>
                        (
                            lMesh.upperAddr(),
                            procFaceMap[proci]
                        )
                    << procFaceCells[proci];
            }
        }

        {
            UIPstream is(Pstream::master(), pBufs);
            label nCells(readLabel(is));
            labelList lowerAddr(is);
            labelList upperAddr(is);
            labelListList faceCells(is);

            PtrList<const lduInterface> interfaces(faceCells.size());
            label inti = 0;
            forAll(faceCells, i)
            {
                if (faceCells[i].size())
                {
                    interfaces.set
                    (
                        inti++,
                        new lduPrimitiveInterface
                        (
                            faceCells[i],
                            labelList(0),
                            labelList(0)
                        )
                    );
                }
            }
            interfaces.setSize(inti);

            const lduPrimitiveMesh lMesh
            (
                nCells,
                lowerAddr,
                upperAddr,
                interfaces,
                lduSchedule(0),
                Pstream::worldComm
            );

            //scalarField diag(is);
            //scalarField lower(is);
            //scalarField upper(is);

        }
    }
//
//            lduMatrix lMat(lMesh);
//            lMat.diag() = diag;
//            lMat.lower() = lower;
//            lMat.upper() = upper;
//
//                const lduInterfaceFieldPtrsList scalarInterfaces(0);
//    const FieldField<Field, scalar> interfaceBouCoeffs(0);
//    const FieldField<Field, scalar> interfaceIntCoeffs(0);
//
//
//
//        }
//
//
//
//        if (Pstream::master())
//        {




    return 0;
}


// ************************************************************************* //
