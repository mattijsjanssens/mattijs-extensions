/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Test-processorAgglomeration

Description
- try and move some bits of lduMesh across and re-stitch
- only coincident cyclics and processor interfaces are stitchable
- all other interfaces get preserved

\*---------------------------------------------------------------------------*/

#include "fvcCellReduce.H"
#include "lduPrimitiveMesh.H"
#include "lduPrimitiveInterface.H"
#include "fvCFD.H"
#include "processorGAMGInterface.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeMatrix(const fvScalarMatrix& mat)
{
    const volScalarField& psi = mat.psi();
    const fvMesh& mesh = psi.mesh();

    {
        Info<< "Writing source" << endl;
        volScalarField source
        (
            IOobject
            (
                "source",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", mat.dimensions(), 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        source.ref().Field<scalar>::operator=(mat.source());
        source.correctBoundaryConditions();
        source.write();
    }

    if (mat.hasDiag())
    {
        Info<< "Writing diag" << endl;
        volScalarField diag
        (
            IOobject
            (
                "diag",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        diag.ref().Field<scalar>::operator=(mat.diag());
        diag.correctBoundaryConditions();
        diag.write();
    }

    if (mat.hasUpper())
    {
        Info<< "Writing upper" << endl;
        surfaceScalarField upper
        (
            IOobject
            (
                "upper",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );
        upper.ref().Field<scalar>::operator=(mat.upper());
        upper.write();

        volScalarField minUpper
        (
            "minUpper",
            fvc::cellReduce(upper, minEqOp<scalar>(), GREAT)
        );
        minUpper.write();

        volScalarField maxUpper
        (
            "maxUpper",
            fvc::cellReduce(upper, maxEqOp<scalar>(), -GREAT)
        );
        maxUpper.write();
    }

    if (mat.hasLower())
    {
        Info<< "Writing lower" << endl;
        surfaceScalarField lower
        (
            IOobject
            (
                "lower",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );
        lower.ref().Field<scalar>::operator=(mat.lower());
        lower.write();

        volScalarField minLower
        (
            "minLower",
            fvc::cellReduce(lower, minEqOp<scalar>(), GREAT)
        );
        minLower.write();

        volScalarField maxLower
        (
            "maxLower",
            fvc::cellReduce(lower, maxEqOp<scalar>(), -GREAT)
        );
        maxLower.write();
    }
}
void printSubMeshInfo
(
    const primitiveMesh& mesh,
    const labelList& cellMap,
    const labelList& faceMap,
    const lduMesh& subMesh
)
{
    forAll(cellMap, cellI)
    {
        Pout<< "cell:" << cellI << " at:" << mesh.cellCentres()[cellMap[cellI]]
            << endl;
    }
    forAll(faceMap, faceI)
    {
        Pout<< "face:" << faceI << " at:" << mesh.faceCentres()[faceMap[faceI]]
            << " lower:" << subMesh.lduAddr().lowerAddr()[faceI]
            << " upper:" << subMesh.lduAddr().upperAddr()[faceI]
            << endl;

    }

    lduInterfacePtrsList ifs(subMesh.interfaces());

    forAll(ifs, patchI)
    {
        if (ifs.set(patchI))
        {
            Pout<< "    patch:" << patchI
                << " fc:" << ifs[patchI].faceCells() << endl;
        }
    }
}
void sendMesh(const lduMesh& mesh, Ostream& os)
{
    const lduAddressing& addressing = mesh.lduAddr();
    lduInterfacePtrsList interfaces(mesh.interfaces());
    boolList validInterface(interfaces.size());
    forAll(interfaces, intI)
    {
        validInterface[intI] = interfaces.set(intI);
    }

    os  << addressing.size()
        << addressing.lowerAddr()
        << addressing.upperAddr()
        << validInterface;

    forAll(interfaces, intI)
    {
        if (interfaces.set(intI))
        {
            os << interfaces[intI].type();
        }
    }
}
autoPtr<lduPrimitiveMesh> receiveMesh(const label comm, Istream& is)
{
    label nCells = readLabel(is);
    labelList lowerAddr(is);
    labelList upperAddr(is);
    boolList validInterface(is);


    // Construct mesh without interfaces
    autoPtr<lduPrimitiveMesh> meshPtr
    (
        new lduPrimitiveMesh
        (
            nCells,
            lowerAddr,
            upperAddr,
            comm,
            true    // reuse
        )
    );

    // Construct GAMGInterfaces
    lduInterfacePtrsList newInterfaces(validInterface.size());
    forAll(validInterface, intI)
    {
        if (validInterface[intI])
        {
            word coupleType(is);

            Pout<< "Received coupleType:" << coupleType << endl;

            newInterfaces.set
            (
                intI,
                // GAMGInterface::New
                // (
                //     coupleType,
                //     intI,
                //     otherMeshes[i-1].rawInterfaces(),
                //     is
                // ).ptr()
                new lduPrimitiveInterface(labelList(0))
            );
        }
    }

    meshPtr().addInterfaces
    (
        newInterfaces,
        lduPrimitiveMesh::nonBlockingSchedule<lduPrimitiveInterface>
        (
            newInterfaces
        )
    );

    return meshPtr;
}


// Temporary replacement for lduPrimitiveMesh::gather
void gather
(
    const label comm,
    const lduMesh& mesh,
    const labelList& procIDs,
    PtrList<lduPrimitiveMesh>& otherMeshes
)
{
    // Force calculation of schedule (since does parallel comms)
    (void)mesh.lduAddr().patchSchedule();

    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        otherMeshes.setSize(procIDs.size()-1);

        // Slave meshes
        for (label i = 1; i < procIDs.size(); i++)
        {
            //Pout<< "on master :"
            //    << " receiving from slave " << procIDs[i] << endl;

            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,          // bufSize
                Pstream::msgType(),
                comm
            );

            label nCells = readLabel(fromSlave);
            labelList lowerAddr(fromSlave);
            labelList upperAddr(fromSlave);
            boolList validInterface(fromSlave);


            // Construct mesh without interfaces
            otherMeshes.set
            (
                i-1,
                new lduPrimitiveMesh
                (
                    nCells,
                    lowerAddr,
                    upperAddr,
                    comm,
                    true    // reuse
                )
            );

            // Construct GAMGInterfaces
            lduInterfacePtrsList newInterfaces(validInterface.size());
            forAll(validInterface, intI)
            {
                if (validInterface[intI])
                {
                    word coupleType(fromSlave);

                    Pout<< "Received coupleType:" << coupleType << endl;

                    newInterfaces.set
                    (
                        intI,
                        // GAMGInterface::New
                        // (
                        //     coupleType,
                        //     intI,
                        //     otherMeshes[i-1].rawInterfaces(),
                        //     fromSlave
                        // ).ptr()
                        new lduPrimitiveInterface(labelList(0))
                    );
                }
            }

            otherMeshes[i-1].addInterfaces
            (
                newInterfaces,
                lduPrimitiveMesh::nonBlockingSchedule<lduPrimitiveInterface>
                (
                    newInterfaces
                )
            );
       }
    }
    else if (findIndex(procIDs, Pstream::myProcNo(comm)) != -1)
    {
        // Send to master

        const lduAddressing& addressing = mesh.lduAddr();
        lduInterfacePtrsList interfaces(mesh.interfaces());
        boolList validInterface(interfaces.size());
        forAll(interfaces, intI)
        {
            validInterface[intI] = interfaces.set(intI);
        }

        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            Pstream::msgType(),
            comm
        );

        toMaster
            << addressing.size()
            << addressing.lowerAddr()
            << addressing.upperAddr()
            << validInterface;

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                toMaster << interfaces[intI].type();
            }
        }
    }
}


labelList upperTriOrder
(
    const label nCells,
    const labelUList& lower,
    const labelUList& upper
)
{
    labelList nNbrs(nCells, 0);

    // Count number of upper neighbours
    forAll(lower, facei)
    {
        if (upper[facei] < lower[facei])
        {
            FatalErrorInFunction
                << "Problem at face:" << facei
                << " lower:" << lower[facei]
                << " upper:" << upper[facei]
                << exit(FatalError);
        }
        nNbrs[lower[facei]]++;
    }

    // Construct cell-upper cell addressing
    labelList offsets(nCells+1);
    offsets[0] = 0;
    forAll(nNbrs, celli)
    {
        offsets[celli+1] = offsets[celli]+nNbrs[celli];
    }

    nNbrs = offsets;

    labelList cellToFaces(offsets.last());
    forAll(upper, facei)
    {
        label celli = lower[facei];
        cellToFaces[nNbrs[celli]++] = facei;
    }

    // Sort

    labelList oldToNew(lower.size());

    labelList order;
    labelList nbr;

    label newFacei = 0;

    for (label celli = 0; celli < nCells; celli++)
    {
        label startOfCell = offsets[celli];
        label nNbr = offsets[celli+1] - startOfCell;

        nbr.setSize(nNbr);
        order.setSize(nNbr);
        forAll(order, i)
        {
            nbr[i] = upper[cellToFaces[offsets[celli]+i]];
        }
        sortedOrder(nbr, order);

        forAll(order, i)
        {
            label index = order[i];
            oldToNew[cellToFaces[startOfCell + index]] = newFacei++;
        }
    }

    return oldToNew;
}


// (local!) subsetting for lduPrimitiveMesh. Adds exposed faces to added
// lduInterface.
autoPtr<lduPrimitiveMesh> subset
(
    const lduMesh& mesh,
    const labelList& region,
    const label currentRegion,
    labelList& cellMap,
    labelList& faceMap,
    labelList& patchMap
)
{
    const labelUList& lower = mesh.lduAddr().lowerAddr();
    const labelUList& upper = mesh.lduAddr().upperAddr();


    // Assign coarse cell numbers. Assign in increasing cell numbering
    // to keep orientation simple.

    labelList reverseCellMap(mesh.lduAddr().size(), -1);
    cellMap.setSize(mesh.lduAddr().size());

    label coarseCellI = 0;
    forAll(region, cellI)
    {
        if (region[cellI] == currentRegion)
        {
            reverseCellMap[cellI] = coarseCellI;
            cellMap[coarseCellI++] = cellI;
        }
    }

    cellMap.setSize(coarseCellI);


    // Detect exposed faces

    labelList reverseFaceMap(lower.size(), -1);
    faceMap.setSize(lower.size());

    DynamicList<label> exposedFaces(lower.size());
    DynamicList<bool> exposedFaceFlips(lower.size());
    DynamicList<label> exposedCells(lower.size());

    label coarseFaceI = 0;
    forAll(lower, faceI)
    {
        label lRegion = region[lower[faceI]];
        label uRegion = region[upper[faceI]];

        if (lRegion == currentRegion)
        {
            if (uRegion == currentRegion)
            {
                // Keep face
                reverseFaceMap[faceI] = coarseFaceI;
                faceMap[coarseFaceI++] = faceI;
            }
            else
            {
                // Upper deleted
                exposedFaces.append(faceI);
                exposedFaceFlips.append(false);
                exposedCells.append(lower[faceI]);
            }
        }
        else if (uRegion == currentRegion)
        {
            // Lower deleted
            exposedFaces.append(faceI);
            exposedFaceFlips.append(true);
            exposedCells.append(upper[faceI]);
        }
    }

    faceMap.setSize(coarseFaceI);


    labelList subLower(faceMap.size());
    labelList subUpper(faceMap.size());
    forAll(faceMap, subFaceI)
    {
        label faceI = faceMap[subFaceI];
        subLower[subFaceI] = reverseCellMap[lower[faceI]];
        subUpper[subFaceI] = reverseCellMap[upper[faceI]];
    }

    const labelList oldToNewFaces
    (
        upperTriOrder
        (
            cellMap.size(),
            subLower,
            subUpper
        )
    );

    labelList orderedLower(UIndirectList<label>(subLower, oldToNewFaces));
    labelList orderedUpper(UIndirectList<label>(subUpper, oldToNewFaces));


    // Copy the interfaces
    lduInterfacePtrsList ifs(mesh.interfaces());

    PtrList<const lduInterface> primitiveInterfaces(ifs.size()+1);

    patchMap.setSize(primitiveInterfaces.size());
    forAll(ifs, patchI)
    {
        patchMap[patchI] = patchI;
        if (ifs.set(patchI))
        {
            //Pout<< "**** Primitive interface at patch " << patchI
            //    << " type:" << ifs[patchI].type() << endl;
            //primitiveInterfaces.set(ifs[patchI].clone());
        }
    }

    // Add a single interface for the exposed cells
    patchMap.last() = -1;
    primitiveInterfaces.set
    (
        primitiveInterfaces.size()-1,
        new lduPrimitiveInterface
        (
            UIndirectList<label>
            (
                reverseCellMap,
                exposedCells
            )()
        )
    );
    Pout<< "Added exposedFaces patch " << primitiveInterfaces.size()-1
        << " with number of cells:" << exposedCells.size() << endl;

    lduInterfacePtrsList newIfs(primitiveInterfaces.size());
    forAll(primitiveInterfaces, i)
    {
        if (primitiveInterfaces.set(i))
        {
            newIfs.set(i, &primitiveInterfaces[i]);
        }
    }

    lduSchedule ps
    (
        lduPrimitiveMesh::nonBlockingSchedule<processorGAMGInterface>(newIfs)
    );

    return autoPtr<lduPrimitiveMesh>
    (
        new lduPrimitiveMesh
        (
            cellMap.size(),
            orderedLower,
            orderedUpper,
            primitiveInterfaces,
            ps,
            Pstream::worldComm
        )
    );
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    //#include "createFields.H"
    //fvScalarMatrix Teqn(fvm::laplacian(DT, T));
    //writeMatrix(Teqn);


    labelList destination(mesh.nCells());

    if (Pstream::master())
    {
        destination = 1;
        for (label cellI = 0; cellI < mesh.nCells()/2; cellI++)
        {
            destination[cellI] = 0;
        }
        Pout<< "Moving " << findIndices(destination, 1).size()
            << " cells to 1" << endl;
    }
    else
    {
        // destination = 0;
        // for (label cellI = 0; cellI < mesh.nCells()/2; cellI++)
        // {
        //     destination[cellI] = 1;
        // }
        destination = 1;

        Pout<< "Moving " << findIndices(destination, 0).size()
            << " cells to 0" << endl;
    }


    const lduMesh& lMesh = mesh;


//XXXXXXXXXXX

// if (false)
// {
//     // Get mesh to keep
// 
//     labelList myCellMap;
//     labelList myFaceMap;
//     labelList myPatchMap;
//     autoPtr<lduPrimitiveMesh> myMeshPtr
//     (
//         subset
//         (
//             lMesh,
//             destination,
//             Pstream::myProcNo(),    // Region to keep
//             myCellMap,
//             myFaceMap,
//             myPatchMap
//         )
//     );
//     const lduPrimitiveMesh& myMesh = myMeshPtr();
//     Pout<< "myMesh:" << myMesh.info() << endl;
//     printSubMeshInfo(mesh, myCellMap, myFaceMap, myMesh);
// 
// 
//     // Where meshes come from
//     labelList procIDs(Pstream::nProcs(Pstream::worldComm));
//     {
//         procIDs[0] = Pstream::myProcNo();
//         label slotI = 1;
//         forAll(procIDs, procI)
//         {
//             if (procI != procIDs[0])
//             {
//                 procIDs[slotI++] = procI;
//             }
//         }
//     }
// 
//     // Get meshes to send
// 
//     PtrList<lduPrimitiveMesh> subMeshes(procIDs.size()-1);
//     labelListList cellMaps(subMeshes.size());
//     labelListList faceMaps(subMeshes.size());
//     labelListList patchMaps(subMeshes.size());
// 
//     forAll(subMeshes, i)
//     {
//         label destProcI = procIDs[i+1];
//         subMeshes.set
//         (
//             i,
//             subset
//             (
//                 lMesh,
//                 destination,
//                 destProcI,   // Region to keep
//                 cellMaps[i],
//                 faceMaps[i],
//                 patchMaps[i]
//             )
//         );
// 
//         const lduPrimitiveMesh& subMesh = subMeshes[i];
//         Pout<< "For destproc:" << destProcI
//             << " subMesh:" << subMesh.info() << endl;
//         printSubMeshInfo(mesh, cellMaps[i], faceMaps[i], subMesh);
//     }
// 
// 
//     labelList cellOffsets;
//     labelList faceOffsets;
//     labelListList faceMap;
//     labelListList boundaryMap;
//     labelListListList boundaryFaceMap;
// 
//     lduPrimitiveMesh allMesh
//     (
//         Pstream::worldComm,
//         identity(Pstream::nProcs(Pstream::worldComm)),
//         labelList(subMeshes.size(), Pstream::myProcNo()),
//         myMesh,
//         subMeshes,
// 
//         cellOffsets,
//         faceOffsets,
//         faceMap,
//         boundaryMap,
//         boundaryFaceMap
//     );
//     Pout<< "allMesh:" << allMesh.info() << endl;
// }
//XXXXXXXXXXX


    const label nProcs = Pstream::nProcs(lMesh.comm());

    // Determine local patches
    boolList isGlobalPatch;
    {
        const lduInterfacePtrsList ifs(lMesh.interfaces());

        isGlobalPatch(ifs.size(), true);
        forAll(ifs, patchi)
        {
            if (ifs[patchi].type() == processorFvPatch::typeName)
            {
                isGlobalPatch[patchi] = false;
            }
        }
    }



    //PtrList<lduPrimitiveMesh> subMeshes(Pstream::nProcs(Pstream::worldComm));
    labelListList cellMaps(nProcs);
    labelListList faceMaps(nProcs);
    labelListList patchMaps(nProcs);


    // Construct local mesh that stays
    autoPtr<lduPrimitiveMesh> myMeshPtr
    (
        subset
        (
            lMesh,
            destination,
            Pstream::myProcNo(),    // Region to keep
            cellMaps[Pstream::myProcNo()],
            faceMaps[Pstream::myProcNo()],
            patchMaps[Pstream::myProcNo()]
        )
    );
    const lduPrimitiveMesh& myMesh = myMeshPtr();


    {
        string oldPrefix = Pout.prefix();
        Pout.prefix() += "myMesh: ";

        Pout<< myMesh.info() << endl;
        printSubMeshInfo
        (
            mesh,
            cellMaps[Pstream::myProcNo()],
            faceMaps[Pstream::myProcNo()],
            myMesh
        );
        Pout.prefix() = oldPrefix;
        Pout<< endl;
    }


    PstreamBuffers pBufs
    (
        Pstream::nonBlocking,
        UPstream::msgType(),
        lMesh.comm()
    );

    for (label destProcI = 0; destProcI < nProcs; destProcI++)
    {
        if (destProcI != Pstream::myProcNo())
        {
            labelList& cMap = cellMaps[destProcI];
            labelList& fMap = faceMaps[destProcI];
            labelList& pMap = patchMaps[destProcI];

            autoPtr<lduPrimitiveMesh> subMeshPtr
            (
                subset
                (
                    lMesh,
                    destination,
                    destProcI,      // Region to keep
                    cMap,
                    fMap,
                    pMap
                )
            );

            const lduPrimitiveMesh& subMesh = subMeshPtr();
            //Pout<< "For destproc:" << destProcI
            //    << " subMesh:" << subMesh.info() << endl;
            //printSubMeshInfo(mesh, cMap, fMap, subMesh);

            // Send subMesh
            UOPstream toDest(destProcI, pBufs);
            sendMesh(subMesh, toDest);
        }
    }


    // Mark all sends as being done
    pBufs.finishedSends();

    // Receive and construct meshes
    PtrList<lduPrimitiveMesh> otherMeshes(nProcs-1);
    label otherI = 0;
    for (label procI = 0; procI < nProcs; procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            UIPstream fromDest(procI, pBufs);
            otherMeshes.set(otherI++,  receiveMesh(lMesh.comm(), fromDest));
            const lduMesh& subMesh = otherMeshes[otherI-1];


            {
                string oldPrefix = Pout.prefix();
                string extra("from " + Foam::name(procI) + ": ");
                Pout.prefix() += extra;

                Pout<< subMesh.info() << endl;
                //printSubMeshInfo(mesh, cMap, fMap, subMesh);

                Pout.prefix() = oldPrefix;
                Pout<< endl;
            }
        }
    }



    Info<< "end" << endl;
}


// ************************************************************************* //
