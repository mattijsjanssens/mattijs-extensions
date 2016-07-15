/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 M Janssens
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
    meshStructure

Description
    Outputs the ordering of the mesh (if any)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "meshStructure.H"
#include "uindirectPrimitivePatch.H"
#include "topoDistanceData.H"
#include "OppositeFaceCellWave.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Determines mesh structure");
    argList::validArgs.append("patches");

    #include "addRegionOption.H"
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    // Get times list
    instantList Times = runTime.times();

    // Set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    #include "createNamedMesh.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Find set of patches from the list of regular expressions provided
    const wordReList patches((IStringStream(args[1])()));
    const labelHashSet patchSet(pbm.patchSet(patches));

    label nFaces = 0;
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        nFaces += pbm[iter.key()].size();
    }

    labelList meshFaces(nFaces);
    nFaces = 0;
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        const polyPatch& pp = pbm[iter.key()];
        forAll(pp, i)
        {
            meshFaces[nFaces++] = pp.start()+i;
        }
    }

    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(mesh.faces(), meshFaces),
        mesh.points()
    );

    //meshStructure ms(mesh, upp);

    List<topoDistanceData> allCellInfo(mesh.nCells());
    List<topoDistanceData> allFaceInfo(mesh.nFaces());

    DynamicList<label> patchFaces(pp.size());
    DynamicList<topoDistanceData> patchData(pp.size());
    forAll(pp, patchFacei)
    {
        patchFaces.append(pp.addressing()[patchFacei]);
        patchData.append(topoDistanceData(patchFacei, 0));
    }

    OppositeFaceCellWave<topoDistanceData> deltaCalc
    (
        mesh,
        patchFaces,
        patchData,
        allFaceInfo,
        allCellInfo,
        mesh.globalData().nTotalCells()+1
    );

    volScalarField fld
    (
        IOobject
        (
            "distance",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("great", dimless, GREAT),
        zeroGradientFvPatchScalarField::typeName
    );
    forAll(fld, celli)
    {
        if (allCellInfo[celli].valid(deltaCalc.data()))
        {
            fld[celli] = allCellInfo[celli].distance();
        }
    }
    fld.correctBoundaryConditions();
    fld.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
