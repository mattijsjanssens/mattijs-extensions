/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 M. Janssens
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

#include "processorColour.H"
#include "processorLduInterface.H"
#include "processorTopologyNew.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorColour, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::processorColour::colour
(
    const lduMesh& mesh,
    labelList& procColour
)
{
    procColour.resize_nocopy(Pstream::nProcs(mesh.comm()));
    procColour = -1;

    // Re-use processor-topology analysis

    const lduInterfacePtrsList patches = mesh.interfaces();

//    // Filter out the non-processor patches
//    DynamicList<label> procPatchIDs;
//    forAll(patches, patchi)
//    {
//        if (patches.set(patchi))
//        {
//            if (isA<processorLduInterface>(patches[patchi]))
//            {
//                procPatchIDs.append(patchi);
//            }
//        }
//    }
//    lduInterfacePtrsList procPatches(procPatchIDs.size());
//    forAll(procPatches, i)
//    {
//        const label patchi = procPatchIDs[i];
//        const auto& pp = patches[patchi];
//        procPatches.set(i, &pp);
//    }
//
//    const processorTopology pt
//    (
//        processorTopology::New<processorLduInterface, lduInterfacePtrsList>
//        (
//            procPatches,
//            mesh.comm()
//        )
//    );


//XXXXX
    labelListList procNeighbours(Pstream::nProcs(mesh.comm()));
    // Fill my entry
    {
        auto& procToProcs = procNeighbours[Pstream::myProcNo(mesh.comm())];
        label n = 0;
        for (const auto& pp : patches)
        {
            if (isA<processorLduInterface>(pp))
            {
                n++;
            }
        }
//        forAll(patches, patchi)
//        {
//            if (patches.set(patchi))
//            {
//                if (isA<processorLduInterface>(patches[patchi]))
//                {
//                    n++;
//                }
//            }
//        }
        procToProcs.setSize(n);
        n = 0;
        for (const auto& pp : patches)
        {
            const auto* ppPtr = isA<processorLduInterface>(pp);
            if (ppPtr)
            {
                procToProcs[n++] = ppPtr->neighbProcNo();
            }
        }
    }
    // Send to master
    Pstream::gatherList(procNeighbours, UPstream::msgType(), mesh.comm());

    DebugVar(procNeighbours);
//XXXXX


    // Use greedy algorithm for now
    // (see e.g. https://iq.opengenus.org/graph-colouring-greedy-algorithm/)

    if (Pstream::master(mesh.comm()))
    {
        DynamicList<label> front;
        front.append(0);    // start from processor 0

        DynamicList<label> newFront;
        while (front.size())
        {
            //Pout<< "Front:" << front << endl;

            newFront.clear();
            for (const label proci : front)
            {
                if (procColour[proci] == -1)
                {
                    //const labelList& nbrs = pt.procNeighbours()[proci];
                    const labelList& nbrs = procNeighbours[proci];
                    const UIndirectList<label> nbrColour(procColour, nbrs);

                    for
                    (
                        label colouri = 0;
                        colouri < Pstream::nProcs();
                        colouri++
                    )
                    {
                        if (!nbrColour.found(colouri))
                        {
                            //Pout<< "Processor:" << proci
                            //    << " allocated colour:" << colouri
                            //    << endl;

                            procColour[proci] = colouri;
                            for (label nbrProci : nbrs)
                            {
                                if (procColour[nbrProci] == -1)
                                {
                                    newFront.append(nbrProci);
                                }
                            }
                            break;
                        }
                    }
                }
            }

            front = std::move(newFront);
        }
    }

    Pstream::broadcast(procColour, mesh.comm());


    //if (false)
    //{
    //    volScalarField volColour
    //    (
    //        IOobject
    //        (
    //            "colour",
    //            mesh.time().timeName(),
    //            mesh,
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE,
    //            false
    //        ),
    //        mesh,
    //        dimensionedScalar(dimless, procColour[Pstream::myProcNo()]),
    //        zeroGradientFvPatchScalarField::typeName
    //    );
    //    volColour.write();
    //}

    const label nColours = max(procColour)+1;

    return nColours;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorColour::processorColour(const lduMesh& mesh)
:
    MeshObject<lduMesh, Foam::MoveableMeshObject, processorColour>(mesh)
{
    nColours_ = colour(mesh, *this);

    Info<< typeName << " : coloured " << Pstream::nProcs()
        << " processors with in total " << nColours_ << " colours" << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

const Foam::processorColour& Foam::processorColour::New(const lduMesh& mesh)
{
    const processorColour* ptr =
        mesh.thisDb().objectRegistry::template cfindObject<processorColour>
        (
            processorColour::typeName
        );

    if (ptr)
    {
        return *ptr;
    }

    processorColour* objectPtr = new processorColour(mesh);

    //regIOobject::store(static_cast<MoveableMeshObject<lduMesh>*>(objectPtr));
    regIOobject::store(objectPtr);

    return *objectPtr;
}


// ************************************************************************* //
