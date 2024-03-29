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

#include "processorTopologyNew.H"
#include "argList.H"
#include "fvMesh.H"
//#include "volFields.H"
#include "processorLduInterface.H"
#include "columnFvMesh.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "processorColour.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//const processorColour& lookup(const lduMesh& mesh)
//{
//    const processorColour* ptr =
//        mesh.thisDb().objectRegistry::template cfindObject<processorColour>
//        (
//            processorColour::typeName
//        );
//
//    if (ptr)
//    {
//        return *ptr;
//    }
//
//    processorColour* objectPtr = new processorColour(mesh);
//
//    regIOobject::store(static_cast<MoveableMeshObject<lduMesh>*>(objectPtr));
//
//    return *objectPtr;
//}


// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const lduMesh& lm = mesh;


    const auto& pc = processorColour::New(lm);
    //const auto& pc = lookup(lm);

    Pout<< "my colour:" << pc.myColour() << endl;
    Pout<< "all colours:" << pc.nColours();

    return 0;



    const lduInterfacePtrsList patches = lm.interfaces();

    // Filter out the non-processor patches
    DynamicList<label> procPatchIDs;
    forAll(patches, patchi)
    {
        if (patches.set(patchi))
        {
            if (isA<processorLduInterface>(patches[patchi]))
            {
                procPatchIDs.append(patchi);
            }
        }
    }
    lduInterfacePtrsList procPatches(procPatchIDs.size());
    forAll(procPatches, i)
    {
        const label patchi = procPatchIDs[i];
        const auto& pp = patches[patchi];
        procPatches.set(i, &pp);
    }

    const processorTopology pt
    (
        processorTopology::New<processorLduInterface, lduInterfacePtrsList>
        (
            procPatches,
            UPstream::worldComm
        )
    );


    Pout<< "procNeighbours:"
        << flatOutput(pt.procNeighbours()[Pstream::myProcNo()])
        << endl;


    labelList colour(Pstream::nProcs(), -1);

    if (Pstream::master())
    {
        DynamicList<label> front;
        front.append(0);    // start from processor 0

        DynamicList<label> newFront;
        while (front.size())
        {
            Pout<< "Front:" << front << endl;

            newFront.clear();
            for (const label proci : front)
            {
                if (colour[proci] == -1)
                {
                    const labelList& nbrs = pt.procNeighbours()[proci];
                    const UIndirectList<label> nbrColour(colour, nbrs);

                    for
                    (
                        label colouri = 0;
                        colouri < Pstream::nProcs();
                        colouri++
                    )
                    {
                        if (!nbrColour.found(colouri))
                        {
                            Pout<< "Processor:" << proci
                                << " allocated colour:" << colouri
                                << endl;

                            colour[proci] = colouri;
                            for (label nbrProci : nbrs)
                            {
                                if (colour[nbrProci] == -1)
                                {
                                    newFront.append(nbrProci);
                                }
                            }
                            break;
                        }
                    }
                }
            }

            Pout<< "    newFront:" << flatOutput(newFront) << endl;

            front = std::move(newFront);
        }
    }

    Pstream::broadcast(colour);


    {
        volScalarField volColour
        (
            IOobject
            (
                "colour",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, colour[Pstream::myProcNo()]),
            zeroGradientFvPatchScalarField::typeName
        );
        //forAll(colour, celli)
        //{
        //    volColour[celli] = colour[celli];
        //}
        //volColour.correctBoundaryConditions();
        volColour.write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
