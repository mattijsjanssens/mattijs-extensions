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

#include "lduInterfaceField.H"
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
#include "lduPrimitiveProcessorInterfaceField.H"

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
    PtrList<lduInterface>& interfaces,

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


class scalarFieldINew
{
public:

    //- Construct null
    scalarFieldINew()
    {}

    //- Construct from Istream
    autoPtr<scalarField> operator()(Istream& is) const
    {
        return autoPtr<scalarField>(new scalarField(is));
    }
};


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

    // ldu mesh
    autoPtr<lduPrimitiveMesh> masterMeshPtr;
    // complete matrix and source
    autoPtr<scalarField> masterSourcePtr;
    autoPtr<lduMatrix> masterMatrixPtr;
    PtrList<lduInterfaceField> masterInterfaces;
    //lduInterfaceFieldPtrsList masterInterfaces;
    FieldField<Field, scalar> masterInterfaceBouCoeffs;
    FieldField<Field, scalar> masterInterfaceIntCoeffs;


    // From processor+local cell to global cell
    labelListList procCellMap(Pstream::nProcs());
    // From my internal face to global internal face
    labelListList procFaceMap(Pstream::nProcs());

    // Per interface face original index in interface
    labelListListList procPatchFaceMap(Pstream::nProcs());


    if (Pstream::master())
    {
        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;

        // Addressing
        label nCells;
        labelList lowerAddr;
        labelList upperAddr;
        PtrList<lduInterface> interfaces;

        // Coefficients
        scalarField diag;
        scalarField lower;
        scalarField upper;
        scalarField source;

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
            masterSourcePtr.reset(new scalarField("value", sourceDict, nCells));
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
        masterMatrixPtr.reset(new lduMatrix(mesh));
        lduMatrix& masterMatrix = masterMatrixPtr();
        masterMatrix.diag().transfer(diag);
        masterMatrix.lower().transfer(lower);
        masterMatrix.upper().transfer(upper);


        // No boundary handling
        masterInterfaces.clear();
        masterInterfaceBouCoeffs.clear();
        masterInterfaceIntCoeffs.clear();


        // Matrix
        // ~~~~~~
        //
        // At this point we have a mesh (lduPrimitiveMesh), matrix
        // coefficients (lduMatrix) and source.
        //


        // Solve
        {
            lduInterfaceFieldPtrsList ifs(masterInterfaces.size());
            forAll(masterInterfaces, i)
            {
                if (masterInterfaces.set(i))
                {
                    ifs.set(i, &masterInterfaces[i]);
                }
            }

            scalarField psi(nCells, 123.0);

            solverPerformance solverPerf = lduMatrix::solver::New
            (
                "p",
                masterMatrix,
                masterInterfaceBouCoeffs,
                masterInterfaceIntCoeffs,
                ifs,    //masterInterfaces,
                solverControls
            )->solve(psi, masterSourcePtr());

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
            List<DynamicList<label>> dynProcCellMap(Pstream::nProcs());
            // Local cell number
            labelList localCell(decomposition.size());
            forAll(decomposition, celli)
            {
                localCell[celli] = dynProcCellMap[decomposition[celli]].size();
                dynProcCellMap[decomposition[celli]].append(celli);
            }
        
            // From my processor to global face
            List<DynamicList<label>> dynProcFaceMap(Pstream::nProcs());
        
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
                    dynProcFaceMap[ownProc].append(facei);
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

            procCellMap.setSize(dynProcCellMap.size());
            procFaceMap.setSize(dynProcFaceMap.size());
            forAll(dynProcCellMap, proci)
            {
                procCellMap[proci].transfer(dynProcCellMap[proci]);
                procFaceMap[proci].transfer(dynProcFaceMap[proci]);
            }

        
        
            // Send subset of mesh
            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
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

                PtrList<lduInterface> interfaces(validNbrs.size());
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
                            tensorField(0),     //tensorField(1, tensor::I),
                            Pstream::msgType()
                        )
                    );
                }
        
                lduPrimitiveMesh procMesh
                (
                    nCells,
                    subLower,
                    subUpper,
                    Pstream::worldComm,
                    true
                );
                procMesh.addInterfaces(interfaces);

                UOPstream os(proci, pBufs);
                lduPrimitiveMeshTools::writeData(procMesh, os);
            }
        }
    }


    // Decomposed data
    // ~~~~~~~~~~~~~~~

    // ldu mesh
    autoPtr<lduPrimitiveMesh> localMeshPtr;
    // complete matrix and source
    autoPtr<scalarField> localSourcePtr;
    autoPtr<lduMatrix> localMatrixPtr;
    PtrList<lduInterfaceField> localInterfaces;
    FieldField<Field, scalar> localInterfaceBouCoeffs;
    FieldField<Field, scalar> localInterfaceIntCoeffs;

    // Receive mesh
    if (Pstream::parRun())
    {
        pBufs.finishedSends();
    
        UIPstream is(Pstream::masterNo(), pBufs);
        localMeshPtr= lduPrimitiveMeshTools::readData(is);
        lduMesh::debug = 1;
        Pout<< localMeshPtr().info() << endl;
    }




    // Attempt 2: use the lduPrimitiveMesh functionality
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (false)  //(Pstream::parRun())
    {
        // Cell to original cell
        labelListList procCellMap(Pstream::nProcs());
        // Internal face to original face
        labelListList procFaceMap(Pstream::nProcs());
        // Interface back to original interface
        labelListList procPatchMap(Pstream::nProcs());
        // Per interface face original index in interface
        //labelListListList procPatchFaceMap(Pstream::nProcs());
        // List of exposed internal faces
        //labelListList procExposedFaceMap(Pstream::nProcs());

        pBufs.clear();
        if (Pstream::master())
        {
            const lduPrimitiveMesh& mesh = masterMeshPtr();
            const label nCells = mesh.lduAddr().size();

            labelList offsets(Pstream::nProcs()+1, nCells);
            offsets[0] = 0;
            const globalIndex globalCells(offsets.xfer());

            labelList decomposition(decompose(nCells, Pstream::nProcs()));

            const lduInterfacePtrsList ifs(mesh.interfaces());
            // Assume all interfaces are global ones ...
            boolList isGlobalInterface(ifs.size(), true);

            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                labelList exposedFaces;
                labelList exposedCells;
                autoPtr<lduPrimitiveMesh> procMeshPtr
                (
                    lduPrimitiveMeshTools::subset
                    (
                        mesh.comm(),
                        mesh,
                        globalCells,
                        ifs,
                        isGlobalInterface,

                        decomposition,
                        proci,

                        procCellMap[proci],
                        procFaceMap[proci],
                        procPatchMap[proci],
                        procPatchFaceMap[proci],

                        exposedFaces, //procExposedFaceMap[proci],
                        exposedCells
                    )
                );
                lduPrimitiveMesh& procMesh = procMeshPtr();

                // Add the exposed internal faces as an interface

                labelList exposedGlobalCells(exposedFaces.size());
                labelList exposedGlobalNbrCells(exposedFaces.size());

                const lduAddressing& addr = mesh.lduAddr();
                forAll(exposedFaces, i)
                {
                    label facei = exposedFaces[i];
                    label own = addr.lowerAddr()[facei];
                    label nei = addr.upperAddr()[facei];
                    label procCelli = exposedCells[i];

                    if (procCelli == own)
                    {
                        exposedGlobalCells[i] = globalCells.toGlobal(own);
                        exposedGlobalNbrCells[i] = globalCells.toGlobal(nei);
                    }
                    else if (procCelli == nei)
                    {
                        exposedGlobalCells[i] = globalCells.toGlobal(nei);
                        exposedGlobalNbrCells[i] = globalCells.toGlobal(own);
                    }
                    else
                    {
                        FatalErrorInFunction << "problem" << exit(FatalError);
                    }
                }

                // Append exposed faces to list of interfaces
                PtrList<lduInterface>& procInterfaces =
                    const_cast<PtrList<lduInterface>&>
                    (
                        procMesh.primitiveInterfaces()
                    );
                PtrList<lduInterface> newInterfaces(procInterfaces.size()+1);
                forAll(procInterfaces, inti)
                {
                    if (procInterfaces.set(inti))
                    {
                        newInterfaces.set
                        (
                            inti,
                            procInterfaces.set(inti, nullptr)
                        );
                    }
                }
                newInterfaces.set
                (
                    procInterfaces.size(),
                    new lduPrimitiveInterface
                    (
                        exposedCells,
                        exposedGlobalCells,
                        exposedGlobalNbrCells
                    )
                );
                procMesh.addInterfaces(newInterfaces);


                // Stream mesh

                UOPstream os(proci, pBufs);
                lduPrimitiveMeshTools::writeData(procMesh, os);
            }
        }

        pBufs.finishedSends();


        // Receive my mesh
        {
            UIPstream is(Pstream::masterNo(), pBufs);
            localMeshPtr = lduPrimitiveMeshTools::readData(is);

            Pout<< "*** NOW:" << endl;
            lduMesh::debug = 1;
            Pout<< localMeshPtr().info() << endl;
        }
    }


     // Now decompose the matrix
    if (Pstream::parRun())
    {
        const lduMesh& masterMesh = masterMeshPtr();
        const lduMatrix& masterMatrix = masterMatrixPtr();
        const lduInterfacePtrsList masterIfs(masterMesh.interfaces());

        pBufs.clear();
 
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream os(proci, pBufs);

            // Cells for proci
           const labelList& procCells = procCellMap[proci];

            // Internal faces for proci
            const labelList& procFaces = procFaceMap[proci];

            os  << UIndirectList<scalar>(masterMatrix.diag(), procCells)
                << UIndirectList<scalar>(masterMatrix.upper(), procFaces)
                << UIndirectList<scalar>(masterMatrix.lower(), procFaces);


            boolList validInterfaces(masterInterfaces.size());
            FieldField<Field, scalar> subBouCoeffs(masterInterfaces.size());
            FieldField<Field, scalar> subIntCoeffs(masterInterfaces.size());
            forAll(masterInterfaces, inti)
            {
                if (masterInterfaces.set(inti))
                {
                    validInterfaces[inti] = true;

                    const scalarField& mbc = masterInterfaceBouCoeffs[int];
                    const scalarField& mic = masterInterfaceIntCoeffs[int];

                    //const label origPatch = procPatchMap[proc][inti];
                    // if (origPatch != -1)
                    // {
                    //     // Existing patch
                    //     const labelUList& pfMap =
                    //         procPatchFaceMap[proci][inti];
                    //     subBouCoeffs.set(inti, new scalarField(mbc, pfMap));
                    //     subIntCoeffs.set(inti, new scalarField(mic, pfMap));
                    // }
                    // else
                    {
                        // Processor patch (from internal faces)
                        const lduInterface& pp = masterInterfaces[inti];
                        subBouCoeffs.set(inti, new scalarField(pp.size()));
                        scalarField& bCoeffs = subBouCoeffs[inti];
                        subIntCoeffs.set(inti, new scalarField(pp.size()));
                        scalarField& iCoeffs = subIntCoeffs[inti];

                        const labelList& fMap = exposedFacesMap[proci][inti];
                        forAll(fMap, i)
                        {
                            if (fMap[i] > 0)
                            {
                                label facei = fMap[i]-1;
                                bCoeffs[i] = masterMatrix.upper()[facei];
                                iCoeffs[i] = masterMatrix.lower()[facei];
                            }
                            else
                            {
                                label facei = -fMap[i]-1;
                                bCoeffs[i] = masterMatrix.lower()[facei];
                                iCoeffs[i] = masterMatrix.upper()[facei];
                            }
                        }
                    }
                }
            }
            os  << validInterfaces << subBouCoeffs << subIntCoeffs;
        }

        pBufs.finishedSends();

        {
            const lduPrimitiveMesh& localMesh = localMeshPtr();
            const lduInterfacePtrsList localIfs(localMesh.interfaces());

            UIPstream is(Pstream::masterNo(), pBufs);

            scalarField diag(is);
            scalarField lower(is);
            scalarField upper(is);

            // Construct matrix
            localMatrixPtr.reset(new lduMatrix(localMesh));
            localMatrixPtr().diag() = diag;
            localMatrixPtr().lower() = lower;
            localMatrixPtr().upper() = upper;

            boolList validNbrs(is);

            localInterfaces.setSize(validNbrs.size());
            forAll(validNbrs, i)
            {
                if (validNbrs[i])
                {
                    localInterfaces.set
                    (
                        i,
                        new lduPrimitiveProcessorInterfaceField(localIfs[i])
                    );
                }
            }


            PtrList<scalarField> intCoeffs(is, scalarFieldINew());
            localInterfaceIntCoeffs.transfer(intCoeffs);

            PtrList<scalarField> bouCoeffs(is, scalarFieldINew());
            localInterfaceBouCoeffs.transfer(bouCoeffs);
        }



        // Solve
        {
            const lduPrimitiveMesh& localMesh = localMeshPtr();
            const lduMatrix& localMatrix = localMatrixPtr();

            lduInterfaceFieldPtrsList ifs(localInterfaces.size());
            forAll(localInterfaces, i)
            {
                if (localInterfaces.set(i))
                {
                    ifs.set(i, &localInterfaces[i]);
                }
            }

            scalarField psi(localMesh.lduAddr().size(), 123.0);

            solverPerformance solverPerf = lduMatrix::solver::New
            (
                "p",
                localMatrix,
                localInterfaceBouCoeffs,
                localInterfaceIntCoeffs,
                ifs,    //localInterfaces,
                solverControls
            )->solve(psi, localSourcePtr());

            solverPerf.print(Info.masterStream(localMesh.comm()));

            DebugVar(psi);
        }
    }

    return 0;
}


// ************************************************************************* //
