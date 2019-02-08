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
//#include "lduPrimitiveInterface.H"
#include "lduPrimitiveProcessorInterface.H"
#include "globalIndex.H"
#include "lduPrimitiveMeshTools.H"
#include "lduPrimitiveInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

labelList decompose(const label nCells, const label nProcs)
{
    // Decompose in identical parts
    const label nLocal = nCells / nProcs;
    labelList offsets(nProcs+1);
    offsets[0] = 0;
    for (label proci = 1; proci < nProcs; proci++)
    {
        offsets[proci] = offsets[proci-1]+nLocal;
    }
    // Any remainder into last processor
    offsets.last() = nCells;

    labelList decomposition(nCells);
    for (label proci = 0; proci < nProcs; proci++)
    {
        label size = offsets[proci+1]-offsets[proci];
        SubList<label>(decomposition, size, offsets[proci]) = proci;
    }
    return decomposition;
}


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
                FatalErrorInFunction << "Asymmetric addressing not supported."
                    << " Cannot find face between " << celli
                    << " and " << nbri << exit(FatalError);
            }

            lower[facei] = e.second();
        }
    }
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


    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);


    // Set up solver to use
    dictionary solverControls("Test-lduSolver");
    solverControls.add("solver", "PBiCGStab");
    solverControls.add("preconditioner", "DILU");
    solverControls.add("tolerance", 1e-06);
    solverControls.add("relTol", 0);


    autoPtr<lduPrimitiveMesh> masterMeshPtr;

    if (Pstream::master())
    {
        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;

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


        // Read source
        scalarField totalSource;
        {
            const IOdictionary sourceDict
            (
                IOobject
                (
                    sourceName,
                    runTime.system(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
            totalSource = scalarField("value", sourceDict, nCells);
        }

        // Construct addressing
        masterMeshPtr.reset
        (
            new lduPrimitiveMesh
            (
                nCells,
                lowerAddr,
                upperAddr,
                interfaces,
                lduSchedule(0),
                Pstream::worldComm
            )
        );
        const lduPrimitiveMesh& mesh = masterMeshPtr();

        // Construct matrix
        lduMatrix lMat(mesh);
        lMat.diag() = diag;
        lMat.lower() = lower;
        lMat.upper() = upper;


        // No boundary handling
        const lduInterfaceFieldPtrsList scalarInterfaces(0);
        const FieldField<Field, scalar> interfaceBouCoeffs(0);
        const FieldField<Field, scalar> interfaceIntCoeffs(0);

        // Solve
        {
            scalarField psi(nCells, 123.0);

            solverPerformance solverPerf = lduMatrix::solver::New
            (
                "p",
                lMat,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                scalarInterfaces,
                solverControls
            )->solve(psi, totalSource);

            solverPerf.print(Info.masterStream(mesh.comm()));

            DebugVar(psi);
        }

        Pstream::parRun() = oldParRun;


        // Distribute the mesh
        // ~~~~~~~~~~~~~~~~~~~

        if (Pstream::parRun())
        {
            labelList decomposition(decompose(nCells, Pstream::nProcs()));

            // From processor+local cell to global cell
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

            forAll(mesh.lowerAddr(), facei)
            {
                label ownCelli = mesh.lowerAddr()[facei];
                label ownProc = decomposition[ownCelli];
                label ownProcCelli = localCell[ownCelli];
                label nbrCelli = mesh.upperAddr()[facei];
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


            // Send subset of mesh
            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                UOPstream os(proci, pBufs);

                const label nCells(procCellMap[proci].size());

                // Internal faces for proci
                const labelList& procFaces = procFaceMap[proci];

                labelList subLower(procFaces.size());
                labelList subUpper(procFaces.size());
                forAll(procFaces, i)
                {
                    label facei = procFaces[i];
                    subLower[i] = localCell[mesh.lowerAddr()[facei]];
                    subUpper[i] = localCell[mesh.upperAddr()[facei]];
                }

                // Construct the interface only for its IO
                const List<DynamicList<label>>& pFaceCells =
                    procFaceCells[proci];

                DynamicList<label> validNbrs(Pstream::nProcs());
                forAll(pFaceCells, nbrProci)
                {
                    if (pFaceCells[nbrProci].size())
                    {
                        validNbrs.append(nbrProci);
                    }
                }
                PtrList<lduPrimitiveProcessorInterface> interfaces
                (
                    validNbrs.size()
                );
                forAll(validNbrs, i)
                {
                    label nbrProci = validNbrs[i];
                    interfaces.set
                    (
                        i,
                        new lduPrimitiveProcessorInterface
                        (
                            pFaceCells[nbrProci],
                            mesh.comm(),
                            proci,
                            nbrProci,
                            tensorField(1, tensor::I),
                            Pstream::msgType()
                        )
                    );
                }

Pout<< "Sending to " << proci
    << " nCells:" << nCells
    << " lower:" << subLower
    << " upper:" << subUpper
    << " validNbrs:" << validNbrs
    << endl;

                os  << nCells
                    << subLower
                    << subUpper
                    << validNbrs;
                forAll(interfaces, i)
                {
                    interfaces[i].write(os);
                }
            }
        }
    }

    autoPtr<lduPrimitiveMesh> localMeshPtr;
    if (Pstream::parRun())
    {
        pBufs.finishedSends();

        UIPstream is(Pstream::masterNo(), pBufs);

        const label nCells(readLabel(is));
        labelList lowerAddr(is);
        labelList upperAddr(is);
        labelList validNbrs(is);

        PtrList<const lduInterface> interfaces(validNbrs.size());
        forAll(validNbrs, i)
        {
            interfaces.set
            (
                i,
                new lduPrimitiveProcessorInterface(is)
            );
        }

        localMeshPtr.reset
        (
            new lduPrimitiveMesh
            (
                nCells,
                lowerAddr,
                upperAddr,
                interfaces,
                lduSchedule(0),
                Pstream::worldComm
            )
        );
        lduMesh::debug = 1;
        Pout<< localMeshPtr().info() << endl;
    }




    // Attempt 2: use the lduPrimitiveMesh functionality
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        pBufs.clear();
        if (Pstream::master())
        {
            const lduPrimitiveMesh& mesh = masterMeshPtr();
            const label nCells = mesh.lduAddr().size();

            labelList offsets(Pstream::nProcs()+1, nCells);
            const globalIndex globalCells(offsets.xfer());

            labelList decomposition(decompose(nCells, Pstream::nProcs()));

            const lduInterfacePtrsList ifs(masterMeshPtr().interfaces());
            // Assume all interfaces are global ones ...
            boolList isGlobalInterface(ifs.size(), true);

            PtrList<labelField> interfaceNbrCells(ifs.size());

            labelListList procCellMap(Pstream::nProcs());
            labelListList procFaceMap(Pstream::nProcs());
            labelListList procPatchMap(Pstream::nProcs());

            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                lduPrimitiveMesh procMesh
                (
                    masterMeshPtr().comm(),
                    masterMeshPtr(),
                    globalCells,
                    ifs,
                    isGlobalInterface,
                    interfaceNbrCells,

                    decomposition,
                    proci,

                    procCellMap[proci],
                    procFaceMap[proci],
                    procPatchMap[proci]
                );

                // Adapt the exposed faces to be processor interfaces
                const lduInterface& couples = procMesh.interfaces().last();
                const lduPrimitiveInterface& lpi =
                    refCast<const lduPrimitiveInterface>(couples);

                List<DynamicList<label>> 


                UOPstream os(proci, pBufs);
                lduPrimitiveMeshTools::writeData(procMesh, os);
            }
        }

        pBufs.finishedSends();

        UIPstream is(Pstream::masterNo(), pBufs);
        localMeshPtr = lduPrimitiveMeshTools::readData(is);


        Pout<< "*** NOW:" << endl;
        lduMesh::debug = 1;
        Pout<< localMeshPtr().info() << endl;
    }


//     // Now send over the matrix
//     if (Pstream::parRun())
//     {
//         PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
// 
//         for (label proci = 0; proci < Pstream::nProcs(); proci++)
//         {
//             UOPstream os(proci, pBufs);
// 
//             // Cells for proci
//             const labelList& procCells = procCellMap[proci];
// 
//             // Internal faces for proci
//             const labelList& procFaces = procFaceMap[proci];
// 
// //             labelList subLower(fMap.size());
// //             labelList subUpper(fMap.size());
// //             forAll(fMap, i)
// //             {
// //                 label facei = fMap[i];
// //                 subLower[i] = localCell[lMesh.lowerAddr()[facei]];
// //                 subUpper[i] = localCell[lMesh.upperAddr()[facei]];
// //             }
// 
//             // Construct the interface only for its IO
//             const List<DynamicList<label>>& pFaceCells =
//                 procFaceCells[proci];
// 
//             DynamicList<label> validNbrs(Pstream::nProcs());
//             forAll(pFaceCells, nbrProci)
//             {
//                 if (pFaceCells[nbrProci].size())
//                 {
//                     validNbrs.append(nbrProci);
//                 }
//             }
//             PtrList<lduPrimitiveProcessorInterface> interfaces
//             (
//                 validNbrs.size()
//             );
//             forAll(validNbrs, i)
//             {
//                 label nbrProci = validNbrs[i];
//                 interfaces.set
//                 (
//                     i,
//                     new lduPrimitiveProcessorInterface
//                     (
//                         pFaceCells[nbrProci],
//                         lMesh.comm(),
//                         proci,
//                         nbrProci,
//                         tensorField(1, tensor::I),
//                         Pstream::msgType()
//                     )
//                 );
//             }
// 
// 
// Pout<< "Sending to " << proci
// << " nCells:" << nCells
// << " lower:" << subLower
// << " upper:" << subUpper
// << " validNbrs:" << validNbrs
// << endl;
// 
// 
// 
//             os  << subLower
//                 << subUpper
//                 << validNbrs;
//             forAll(interfaces, i)
//             {
//                 interfaces[i].write(os);
//             }
//         }
// 
//         pBufs.finishedSends();
// 
//         for (label proci = 0; proci < Pstream::nProcs(); proci++)
//         {
//             UIPstream is(proci, pBufs);
// 
//             scalarField diag(is);
//             scalarField lower(is);
//             scalarField upper(is);
// 
//             // Construct matrix
//             lduMatrix lMat(lMesh);
//             lMat.diag() = diag;
//             lMat.lower() = lower;
//             lMat.upper() = upper;
// 
//             labelList validNbrs(is);
// 
//             lduInterfaceFieldPtrsList scalarInterfaces(validNbrs.size());
//             forAll(validNbrs, i)
//             {
//                 const lduInterface& intf = lMesh.interfaces()
// 
//                 scalarInterfaces.set
//                 (
//                     i,
//                     new lduPrimitiveProcessorInterfaceField(is)
//                 );
//             }
//         }
//    }


//         // No boundary handling
//         const lduInterfaceFieldPtrsList scalarInterfaces(0);
//         const FieldField<Field, scalar> interfaceBouCoeffs(0);
//         const FieldField<Field, scalar> interfaceIntCoeffs(0);
// 
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
