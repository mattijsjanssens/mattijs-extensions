/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

\*---------------------------------------------------------------------------*/

#include "lduPrimitiveMeshTools.H"
#include "lduPrimitiveInterface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::lduPrimitiveMeshTools::writeData(const lduMesh& mesh, Ostream& os)
{
    const lduAddressing& addressing = mesh.lduAddr();
    lduInterfacePtrsList interfaces(mesh.interfaces());
    boolList validInterface(interfaces.size());
    forAll(interfaces, intI)
    {
        validInterface[intI] = interfaces.set(intI);
    }

    os
        << mesh.comm()
        << addressing.size()
        << addressing.lowerAddr()
        << addressing.upperAddr()
        << validInterface;

    // Write raw
    forAll(interfaces, intI)
    {
        if (interfaces.set(intI))
        {
            os << interfaces[intI].type();
    //         //refCast<const lduPrimitiveInterface>
    //         //(interfaces[intI]).write(os);
            interfaces[intI].write(os);
        }
    }

    // Write as dictionary
    // {
    //     os  << nl << token::BEGIN_LIST << incrIndent << nl;
    //     forAll(interfaces, intI)
    //     {
    //         if (interfaces.set(intI))
    //         {
    //             os  << indent << token::BEGIN_BLOCK << nl
    //                 << incrIndent;
    //             interfaces[intI].write(os);
    //             os  << decrIndent
    //                 << indent << token::END_BLOCK << endl;
    //         }
    //     }
    //     os  << decrIndent << token::END_LIST;
    // }

    // Check state of IOstream
    os.check("polyBoundaryMesh::writeData(Ostream& os) const");

    return true;
}


Foam::autoPtr<Foam::lduPrimitiveMesh>
Foam::lduPrimitiveMeshTools::readData(Istream& is)
{
    label comm = readLabel(is);
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

    // // Construct GAMGInterfaces
    // lduInterfacePtrsList newInterfaces(validInterface.size());
    // forAll(validInterface, intI)
    // {
    //     if (validInterface[intI])
    //     {
    //         word coupleType(is);
    // 
    //         Pout<< "Received coupleType:" << coupleType << endl;
    // 
    //         newInterfaces.set
    //         (
    //             intI,
    //             // GAMGInterface::New
    //             // (
    //             //     coupleType,
    //             //     intI,
    //             //     otherMeshes[i-1].rawInterfaces(),
    //             //     is
    //             // ).ptr()
    //             new lduPrimitiveInterface(is)
    //         );
    //     }
    // }
    // 
    // meshPtr().addInterfaces
    // (
    //     newInterfaces,
    //     lduPrimitiveMesh::nonBlockingSchedule<lduPrimitiveInterface>
    //     (
    //         newInterfaces
    //     )
    // );

    PtrList<lduInterface> newInterfaces(validInterface.size());
    forAll(validInterface, intI)
    {
        if (validInterface[intI])
        {
            word coupleType(is);
            Pout<< "Received coupleType:" << coupleType << endl;
            newInterfaces.set
            (
                intI,
                lduInterface::New
                (
                    coupleType,
                    is
                )
            );
        }
    }
    meshPtr().addInterfaces(newInterfaces);

    return meshPtr;
}


void Foam::lduPrimitiveMeshTools::swapCellData
(
    const lduInterfacePtrsList& ifs,
    const labelUList& cellData,
    PtrList<labelField>& nbrData
)
{
    nbrData.setSize(ifs.size());

    // Initialise transfer of global cells
    forAll(ifs, patchi)
    {
        if (ifs.set(patchi))
        {
            ifs[patchi].initInternalFieldTransfer
            (
                Pstream::commsTypes::nonBlocking,
                cellData
            );
        }
    }

    if (Pstream::parRun())
    {
        Pstream::waitRequests();
    }

    forAll(ifs, patchi)
    {
        if (ifs.set(patchi))
        {
            nbrData.set
            (
                patchi,
                ifs[patchi].internalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    cellData
                )
            );
        }
    }
}


Foam::autoPtr<Foam::lduPrimitiveMesh>
Foam::lduPrimitiveMeshTools::subset
(
    const label comm,
    const lduMesh& mesh,
    const globalIndex& globalNumbering,
    const lduInterfacePtrsList& interfaces,
    const boolList& isGlobalInterface,

    const labelList& region,
    const label currentRegion,

    labelList& cellMap,
    labelList& faceMap,
    labelList& patchMap,
    labelListList& patchFaceMap,
    labelList& exposedFaceMap,
    labelList& exposedFaceCells
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
    exposedFaceMap.setSize(lower.size());
    exposedFaceCells.setSize(lower.size());
    //DynamicList<label> exposedFaces(lower.size());
    //DynamicList<bool> exposedFaceFlips(lower.size());
    //DynamicList<label> exposedCells(lower.size());
    //DynamicList<label> exposedGlobalCells(lower.size());
    //DynamicList<label> exposedGlobalNbrCells(lower.size());

    label coarseFaceI = 0;
    label exposedFaceI = 0;

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
                exposedFaceMap[exposedFaceI] = faceI;
                exposedFaceCells[exposedFaceI] = reverseCellMap[lower[faceI]];
                exposedFaceI++;
                //exposedFaces.append(faceI);
                //exposedFaceFlips.append(false);
                //exposedCells.append(reverseCellMap[lower[faceI]]);
                //exposedGlobalCells.append
                //(
                //    globalNumbering.toGlobal(lower[faceI])
                //);
                //exposedGlobalNbrCells.append
                //(
                //    globalNumbering.toGlobal(upper[faceI])
                //);
            }
        }
        else if (uRegion == currentRegion)
        {
            // Lower deleted
            exposedFaceMap[exposedFaceI] = faceI;
            exposedFaceCells[exposedFaceI] = reverseCellMap[upper[faceI]];
            exposedFaceI++;
            //exposedFaces.append(faceI);
            //exposedFaceFlips.append(true);
            //exposedCells.append(reverseCellMap[upper[faceI]]);
            //exposedGlobalCells.append
            //(
            //    globalNumbering.toGlobal(upper[faceI])
            //);
            //exposedGlobalNbrCells.append
            //(
            //    globalNumbering.toGlobal(lower[faceI])
            //);
        }
    }

    faceMap.setSize(coarseFaceI);
    exposedFaceMap.setSize(exposedFaceI);
    exposedFaceCells.setSize(exposedFaceI);


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
        lduPrimitiveMesh::upperTriOrder
        (
            cellMap.size(),
            subLower,
            subUpper
        )
    );

    labelList orderedLower(UIndirectList<label>(subLower, oldToNewFaces));
    labelList orderedUpper(UIndirectList<label>(subUpper, oldToNewFaces));


    // Map all the interfaces
    patchFaceMap.setSize(interfaces.size());
    patchMap.setSize(interfaces.size());
    PtrList<lduInterface> newInterfaces(interfaces.size());
    label newInti = 0;

    forAll(interfaces, inti)
    {
        const lduInterface& intf = interfaces[inti];
        const labelUList& fc = intf.faceCells();

        DynamicList<label> patchFaces(fc.size());
        forAll(fc, i)
        {
            if (region[fc[i]] == currentRegion)
            {
                patchFaces.append(i);
            }
        }

        if (isGlobalInterface[inti] || patchFaces.size())
        {
            newInterfaces.set(newInti, intf.clone(newInti, patchFaces));
            patchFaceMap[newInti].transfer(patchFaces);
            patchMap[newInti] = inti;
            newInti++;
        }
    }
    newInterfaces.setSize(newInti);
    patchFaceMap.setSize(newInti);
    patchMap.setSize(newInti);

    lduInterfacePtrsList newIfs(newInterfaces.size());
    forAll(newInterfaces, i)
    {
        if (newInterfaces.set(i))
        {
            newIfs.set(i, &newInterfaces[i]);
        }
    }

    lduSchedule ps
    (
        lduPrimitiveMesh::nonBlockingSchedule<lduPrimitiveInterface>(newIfs)
    );

    return autoPtr<lduPrimitiveMesh>
    (
        new lduPrimitiveMesh
        (
            cellMap.size(),
            orderedLower,
            orderedUpper,
            newInterfaces,
            ps,
            comm
        )
    );
}


// ************************************************************************* //
