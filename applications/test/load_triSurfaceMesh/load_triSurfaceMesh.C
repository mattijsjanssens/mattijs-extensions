/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    snappyHexMesh

Description
    Automatic split hex mesher. Refines and snaps to surface.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "snappyRefineDriver.H"
#include "snappySnapDriver.H"
#include "snappyLayerDriver.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "decompositionMethod.H"
#include "noDecomp.H"
#include "fvMeshDistribute.H"
#include "wallPolyPatch.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "layerParameters.H"
#include "vtkSetWriter.H"
#include "faceSet.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "cellModeller.H"
#include "uindirectPrimitivePatch.H"
#include "surfZoneIdentifierList.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurface.H"
#include "globalIndex.H"
#include "IOmanip.H"
#include "fvMeshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    // Read meshing dictionary
    const word dictName("snappyHexMeshDict");
    #include "setSystemMeshDictionaryIO.H"
    const IOdictionary meshDict(dictIO);

    // all surface geometry
    const dictionary& geometryDict = meshDict.subDict("geometry");

    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            //mesh.time().findInstance("triSurface", word::null),// instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        true    //meshDict.lookupOrDefault("singleRegionName", true)
    );


    forAll(allGeometry, geomi)
    {
        Pout<< "geom:" << allGeometry[geomi].name()
            << " size:" << allGeometry[geomi].size()
            << endl;
    }


    runTime++;
    forAll(allGeometry, geomi)
    {
        allGeometry[geomi].write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
