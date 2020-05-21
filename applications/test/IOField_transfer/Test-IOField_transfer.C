/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 M. Janssens
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
    sphericalTensorFieldTest

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvCFD.H"
//#include "skewCorrectionVectors.H"
//#include "volFields.H"
//#include "surfaceFields.H"
//#include "pisoControl.H"
#include "IOField.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const objectRegistry& patchData
(
    objectRegistry& db,
    const word& world,
    const word& region,
    const word& patch
)
{
    const objectRegistry& worldDb = db.subRegistry
    (
        world,
        true,   //const bool forceCreate = false,
        false   //const bool recursive = false
    );
    const objectRegistry& regionDb = worldDb.subRegistry
    (
        region,
        true,   //const bool forceCreate = false,
        false   //const bool recursive = false
    );
    const objectRegistry& patchDb = regionDb.subRegistry
    (
        patch,
        true,   //const bool forceCreate = false,
        false   //const bool recursive = false
    );

    return patchDb; //const_cast<objectRegistry&>(patchDb);
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    if (UPstream::allWorlds().size() != 2)
    {
        FatalErrorInFunction << "Only can use two worlds" << exit(FatalError);
    }


//    Info<< "Reading field p\n" << endl;
//    volScalarField vfld
//    (
//        IOobject
//        (
//            "p",
//            runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );
//    vfld.dimensions().reset(dimLength);
//    vfld == mag(mesh.C());


    const word remotePatch("movingWall");
    word remoteWorld;
    for (const word& world : UPstream::allWorlds())
    {
        if (world != UPstream::myWorld())
        {
            remoteWorld = world;
            break;
        }
    }

DebugVar(UPstream::myWorldID());
DebugVar(remotePatch);
DebugVar(remoteWorld);
DebugVar(mesh.name());


    // Get remote world/region/patch database
    const objectRegistry& obj = patchData
    (
        runTime,
        remoteWorld,
        mesh.name(),
        remotePatch
    );

    // Put some data into obj
    if (UPstream::myWorldID() == 0)
    {
        autoPtr<IOField<scalar>> fldPtr
        (
            new IOField<scalar>
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    obj,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                )
            )
        );
        IOField<scalar>& fld = fldPtr();
        fld.setSize(3, 123.0);

        IOField<scalar>::store(fldPtr);
    }

    DebugVar(obj.sortedToc());


    // Send over data to correct world
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);


    }

    return 0;
}


// ************************************************************************* //
