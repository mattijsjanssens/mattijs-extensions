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

    forAll(interfaces, intI)
    {
        if (interfaces.set(intI))
        {
            os << interfaces[intI].type();
            refCast<const lduPrimitiveInterface>(interfaces[intI]).write(os);
        }
    }
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
                new lduPrimitiveInterface(is)
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


// ************************************************************************* //
