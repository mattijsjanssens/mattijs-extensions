/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    flatten2DMesh

Group
    grpMeshGenerationUtilities

Description
    Opposite of extrude2DMesh.
    Note: could also be extrude option in extrudeMesh since
    - reads surface
    - outputs a mesh

Note

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "indirectPrimitivePatch.H"
#include "polyTopoChange.H"
#include "PatchTools.H"
#include "topoSet.H"
#include "EdgeMap.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("patches");

    #include "addOverwriteOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const bool overwrite = args.optionFound("overwrite");

    // Keep only the faces on the selected patches. Create cells out of
    // the faces, faces out of the edges.

    const wordReList patchNames(IStringStream(args.args()[1])());

DebugVar(patchNames);


    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    label nBnd = mesh.nFaces()-mesh.nInternalFaces();

    labelHashSet patchSet(pbm.patchSet(patchNames));
    const labelList patchIds(patchSet.sortedToc());

    DynamicList<label> patchFaces;
    PackedBoolList isPatchFace(nBnd);
    {
        forAll(patchIds, i)
        {
            label patchi = patchIds[i];
            const polyPatch& pp = pbm[patchi];

            forAll(pp, i)
            {
                patchFaces.append(pp.start()+i);
                isPatchFace[pp.start()-mesh.nInternalFaces()+i] = true;
            }
        }
    }

    indirectPrimitivePatch allPatch
    (
        IndirectList<face>(mesh.faces(), patchFaces),
        mesh.points()
    );


    polyTopoChange meshMod(pbm.size(), false);

    // For now add in order. Could be renumbered.
    labelList faceToCell(patchFaces.size());
    forAll(faceToCell, i)
    {
        faceToCell[i] = meshMod.addCell
        (
            -1,     //const label masterPointID,
            -1,     //const label masterEdgeID,
            -1,     //const label masterFaceID,
            -1,     //const label masterCellID,
            -1      //const label zoneID
        );
    }

DebugVar(faceToCell);

    // Add points
    const labelList& allPatchPoints = allPatch.meshPoints();

    forAll(allPatchPoints, ppi)
    {
        label pointi = allPatchPoints[ppi];

        Pout<< "Adding point " << mesh.points()[pointi]
            << " from local point:" << ppi
            << " from mesh point " << pointi << endl;
        meshMod.addPoint
        (
            mesh.points()[pointi],
            pointi, //const label masterPointID,
            -1,     //const label zoneID,
            true    //const bool inCell
        );
    }

    // faces
    const labelList own(PatchTools::edgeOwner(allPatch));
    const edgeList& edges = allPatch.edges();
    const labelListList& edgeFaces = allPatch.edgeFaces();

    face f(2);

    for (label edgei = 0; edgei < allPatch.nInternalEdges(); edgei++)
    {
        const labelList& eFaces = edgeFaces[edgei];

        label nei = eFaces[0];
        if (nei == own[edgei])
        {
            nei = eFaces[1];
        }

        const edge& e = edges[edgei];
        f[0] = e[0];
        f[1] = e[1];

        Pout<< "Adding internal face " << f << " own:" << faceToCell[own[edgei]]
            << " nei:" << faceToCell[nei] << endl;

        meshMod.addFace
        (
            f,          //const face& f,
            faceToCell[own[edgei]], //const label own,
            faceToCell[nei],        //const label nei,
            -1,         //const label masterPointID,
            -1,         //const label masterEdgeID,
            -1,         //const label masterFaceID,
            false,      //const bool flipFaceFlux,
            -1,         //const label patchID,
            -1,         //const label zoneID,
            false       //const bool zoneFlip
        );
    }

    // Get the patch from the patches which are not included

    EdgeMap<label> nbrPatchID(nBnd);
    {
        SubList<face> boundaryFaces(mesh.faces(), nBnd, mesh.nInternalFaces());
        primitivePatch boundaryPatch(boundaryFaces, mesh.points());

        const labelListList& edgeFaces = boundaryPatch.edgeFaces();

        forAll(edgeFaces, edgei)
        {
            const edge& e = boundaryPatch.edges()[edgei];
            Pout<< "edge:" << edgei
                << " verts:" << boundaryPatch.localPoints()[e[0]]
                << boundaryPatch.localPoints()[e[1]] << endl;

            const edge meshE
            (
                boundaryPatch.meshPoints()[e[0]],
                boundaryPatch.meshPoints()[e[1]]
            );

            const labelList& eFaces = edgeFaces[edgei];
            if (eFaces.size() >= 2)
            {
                label patch0 = pbm.patchID()[eFaces[0]];
                label patch1 = pbm.patchID()[eFaces[1]];

                Pout<< "    patch0:" << patch0 << endl;
                Pout<< "    patch1:" << patch1 << endl;
                Pout<< "    eFaces:" << eFaces << endl;

                if (patchSet.found(patch0) && !patchSet.found(patch1))
                {
                    Pout<< "    marking face:"
                        << boundaryPatch.faceCentres()[eFaces[0]]
                        << " with patch " << patch1 << endl;
                    nbrPatchID.insert(meshE, patch1);
                }
                else if (patchSet.found(patch1) && !patchSet.found(patch0))
                {
                    Pout<< "    marking face:"
                        << boundaryPatch.faceCentres()[eFaces[1]]
                        << " with patch " << patch0 << endl;
                    nbrPatchID.insert(meshE, patch0);
                }
            }
            else
            {
                FatalErrorInFunction
                    << "edge:" << boundaryPatch.edges()[edgei]
                    << " eFaces:" << eFaces
                    << exit(FatalError);
            }
        }
    }

DebugVar(nbrPatchID);


    for
    (
        label edgei = allPatch.nInternalEdges();
        edgei < allPatch.nEdges();
        edgei++
    )
    {
        const edge& e = edges[edgei];
        const labelList& eFaces = edgeFaces[edgei];
        if (eFaces.size() > 1)
        {
            WarningInFunction
                << "Boundary edge " << edgei << " with vertices "
                << allPatch.points()[allPatch.meshPoints()[e[0]]]
                << allPatch.points()[allPatch.meshPoints()[e[1]]]
                << " has " << eFaces.size() << " patch faces connected"
                << endl;
        }

        const edge meshE
        (
            allPatch.meshPoints()[e[0]],
            allPatch.meshPoints()[e[1]]
        );


        // We should have a nbrPatchID
        label allFacei = eFaces[0];
        label patchi = nbrPatchID[meshE];

        Pout<< "Boundary edge:" << e
            << " verts:" << allPatch.points()[allPatch.meshPoints()[e[0]]]
            << allPatch.points()[allPatch.meshPoints()[e[1]]]
            << " own:" << faceToCell[allFacei]
            << " on patch:" << patchi << endl;

        f[0] = e[0];
        f[1] = e[1];
        meshMod.addFace
        (
            f,          //const face& f,
            faceToCell[allFacei], //const label own,
            -1,         //const label nei,
            -1,         //const label masterPointID,
            -1,         //const label masterEdgeID,
            -1,         //const label masterFaceID,
            false,      //const bool flipFaceFlux,
            patchi,     //const label patchID,
            -1,         //const label zoneID,
            false       //const bool zoneFlip
        );
    }

    autoPtr<fvMesh> newMeshPtr;
    autoPtr<mapPolyMesh> mapPtr
    (
        meshMod.makeMesh
        (
            newMeshPtr,
            mesh,
            mesh
        )
    );


    if (!overwrite)
    {
        runTime++;
        newMeshPtr().setInstance(runTime.timeName());
    }
    else
    {
        newMeshPtr().setInstance("constant");
    }

    Info<< "\nWriting mesh to time = " << runTime.timeName() << nl << endl;

    newMeshPtr().write();
    topoSet::removeFiles(mesh);
    //processorMeshes::removeFiles(mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
