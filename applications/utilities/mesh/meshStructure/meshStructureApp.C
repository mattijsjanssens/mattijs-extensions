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
//#include "IOobjectList.H"
#include "fvMesh.H"
// #include "polyTopoChange.H"
// #include "ReadFields.H"
#include "volFields.H"
// #include "surfaceFields.H"
// #include "SortableList.H"
// #include "decompositionMethod.H"
// #include "renumberMethod.H"
#include "zeroGradientFvPatchFields.H"
// #include "CuthillMcKeeRenumber.H"
#include "meshStructure.H"
// #include "cellSet.H"
// #include "faceSet.H"
// #include "pointSet.H"
#include "uindirectPrimitivePatch.H"

using namespace Foam;


// Create named field from labelList for postprocessing
tmp<volScalarField> createScalarField
(
    const fvMesh& mesh,
    const word& name,
    const labelList& elems
)
{
    tmp<volScalarField> tfld
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& fld = tfld.ref();

    forAll(fld, celli)
    {
       fld[celli] = elems[celli];
    }

    return tfld;
}


void opposingFaceLabels
(
    const primitiveMesh& mesh,
    const label cellI,
    const label masterFaceLabel,
    DynamicList<label>& oppositeFaceLabels
)
{
    // Variant of cell::opposingFaceLabel

    // Algorithm:
    // Go through all the faces of the cell and find the one which
    // does not share a single vertex with the master face.  If there
    // are two or more such faces, return the first one and issue a
    // warning; if there is no opposite face, return -1;

    const face& masterFace = mesh.faces()[masterFaceLabel];

    const labelList& curFaceLabels = mesh.cells()[cellI];

    oppositeFaceLabels.clear();

    forAll(curFaceLabels, facei)
    {
        // Compare the face with the master
        const face& curFace = mesh.faces()[curFaceLabels[facei]];

        // Skip the master face
        if (curFaceLabels[facei] != masterFaceLabel)
        {
            bool sharedPoint = false;

            // Compare every vertex of the current face agains the
            // vertices of the master face
            forAll(curFace, pointi)
            {
                const label l = curFace[pointi];

                forAll(masterFace, masterPointi)
                {
                    if (masterFace[masterPointi] == l)
                    {
                        sharedPoint = true;
                        break;
                    }
                }

                if (sharedPoint) break;
            }

            // If no points are shared, this is the opposite face
            if (!sharedPoint)
            {
                // Found opposite face
                oppositeFaceLabel.append(curFaceLabels[facei]);
            }
        }
    }
}


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
    List<topoDistanceData> allFaceInfo(mesh.nCells());
    {
        // Start of changes
        DynamicList<label> frontFaces(pp.size());
        DynamicList<topoDistanceData> frontData(pp.size());
        forAll(pp, patchFacei)
        {
            label meshFacei = pp.addressing()[patchFacei];

            // Make sure face present only once in initial front
            if (!allFaceInfo[meshFacei].valid())
            {
                allFaceInfo[meshFacei] = topoDistanceData(patchFacei, 0);
                frontData.append(allFaceInfo[meshFacei]);
                frontFaces.append(meshFacei);
            }
        }


        boolList isNewFrontFace(mesh.nFaces());

        int td;
        DynamicList<label> oppositeFaceLabels;

        while (true)
        {
            DynamicList<label> newFrontFaces(frontFaces.size());
            isNewFrontFace = false;

            // Collect opposite face
            forAll(frontFaces, i)
            {
                label facei = frontFaces[i];

                {
                    // Owner side

                    label own = mesh.faceOwner()[facei];
                    opposingFaceLabels
                    (
                        mesh,
                        own,
                        facei,
                        oppositeFaceLabels
                    );
                    if (oppositeFaceLabels.size())
                    {
                        bool propagate = allCellInfo[own].updateCell
                        (
                            mesh,
                            own,
                            facei,
                            allFaceInfo[facei],
                            1e-6,                   //tol,
                            td
                        );

                        if (propagate && oppositeFaceLabels.size() == 1)
                        {
                            label oppFacei = oppositeFaceLabels[0];

                            bool propagate = allFaceInfo[oppFacei].updateFace
                            (
                                mesh,
                                oppFacei,
                                allCellInfo[own],
                                1e-6,                   //tol,
                                td
                            );
                            if (propagate && isNewFrontFace.set(oppFacei))
                            {
                                newFrontFaces.append(oppFacei);
                            }
                        }
                    }
                }
                if (mesh.isInternalFace(facei))
                {
                    label nei = mesh.faceNeighbour()[facei];
                    opposingFaceLabels
                    (
                        mesh,
                        nei,
                        facei,
                        oppositeFaceLabels
                    );
                    if (oppositeFaceLabels.size())
                    {
                        int td;
                        bool propagate = allCellInfo[nei].updateCell
                        (
                            mesh,
                            nei,
                            facei,
                            allFaceInfo[facei],
                            1e-6,                   //tol,
                            td
                        );

                        if (propagate && oppositeFaceLabels.size() == 1)
                        {
                            label oppFacei = oppositeFaceLabels[0];

                            bool propagate = allFaceInfo[oppFacei].updateFace
                            (
                                mesh,
                                oppFacei,
                                allCellInfo[nei],
                                1e-6,                   //tol,
                                td
                            );
                            if (propagate && isNewFrontFace.set(oppFacei))
                            {
                                newFrontFaces.append(oppFacei);
                            }
                        }
                    }
                }
            }


            //TBD. Synchronise across coupled faces ...

            Pout<< "old front:" << front.size()
                << " new front:" << newFrontFaces.size() << endl;

            if (returnReduce(newFrontFaces.size(), sumOp<label>()) == 0)
            {
                break;
            }

            frontFaces.transfer(newFrontFaces);
        }
    }


//     createScalarField
//     (
//         mesh,
//         "cellLayer",
//         ms.cellLayer()
//     )().write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
