/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    integrateWallDist

Description
    Integrate surfaceScalarField away from patch

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "fvMesh.H"
#include "wallPointData.H"
#include "FaceCellWave.H"
#include "fvcGrad.H"
#include "localMin.H"
#include "localMax.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("patch");

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();
    #include "createMesh.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Find set of patches from the list of regular expressions provided
    const word patchName((IStringStream(args[1])()));
    label patchi = pbm.findPatchID(patchName);

    const fvPatch& fvp = mesh.boundary()[patchi];
    const vectorField& Cf = fvp.Cf();
    const polyPatch& pp = fvp.patch();
    const labelListList& edgeFaces = pp.edgeFaces();
    const labelListList& faceEdges = pp.faceEdges();

    // Addressing: owner = geometrically lower face, neighbour = geometrically
    // higher face
    labelList edgeOwner(pp.nEdges(), -1);
    labelList edgeNeighbour(pp.nEdges(), -1);
    labelList nInEdges(pp.size(), 0);
    {
        forAll(edgeFaces, edgei)
        {
            const edge& e = pp.edges()[edgei];
            const point& p0 = pp.localPoints()[e[0]];
            const point& p1 = pp.localPoints()[e[1]];
            point edgeMid(0.5*(p0+p1));

            const labelList& eFaces = edgeFaces[edgei];
            forAll(eFaces, i)
            {
                label facei = eFaces[i];
                const point& fc = Cf[facei];

                if (fc.x() < edgeMid.x())
                {
                    edgeOwner[edgei] = facei;
                }
                else
                {
                    edgeNeighbour[edgei] = facei;

                    if (eFaces.size() > 1)
                    {
                        nInEdges[facei]++;
                    }
                }
            }
        }
    }

DebugVar(edgeOwner);
DebugVar(edgeNeighbour);
DebugVar(nInEdges);




    // Start seed: lowest edges

    DynamicList<label> frontEdges;
    DynamicList<Pair<scalar>> frontInfo;
    {
        forAll(pp.edges(), edgei)
        {
            const edge& e = pp.edges()[edgei];
            const point& p0 = pp.localPoints()[e[0]];
            const point& p1 = pp.localPoints()[e[1]];

            if (p0.x() < 1e-5 && p1.x() < 1e-5)
            {
                Pout<< "Front edge:"
                    << pp.localPoints()[e[0]]
                    << pp.localPoints()[e[1]]
                    << endl;

                frontEdges.append(edgei);
                frontInfo.append(Pair<scalar>(mag(p0-p1), 0.0));
            }
        }
    }

    List<Pair<scalar>> allFaceInfo(pp.size(), Pair<scalar>(0.0, 0.0));
    List<Pair<scalar>> allEdgeInfo(pp.nEdges(), Pair<scalar>(0.0, 0.0));
    forAll(frontEdges, i)
    {
        allEdgeInfo[frontEdges[i]] = frontInfo[i];
    }

    while (true)
    {
        // New front
        PackedBoolList isNewFrontEdge(pp.nEdges());
        DynamicList<label> newFrontEdges;

DebugVar(frontEdges.size());

        // Check front, accumulate
        forAll(frontEdges, fronti)
        {
            label edgei = frontEdges[fronti];
            const edge& e = pp.edges()[edgei];

            const labelList& eFaces = edgeFaces[edgei];

            forAll(eFaces, eFacei)
            {
                label facei = eFaces[eFacei];
                if (facei == edgeNeighbour[edgei] && nInEdges[facei] > 0)
                {
                    Pout<< "Updating face:" << Cf[facei]
                        << " with info from edge:" << pp.localPoints()[e[0]]
                        << pp.localPoints()[e[1]] << endl;

                    // Accumulate
                    allFaceInfo[facei].first() += allEdgeInfo[edgei].first();
                    allFaceInfo[facei].second() += allEdgeInfo[edgei].second();
                    nInEdges[facei]--;

                    if (nInEdges[facei] == 0)
                    {
                        Pout<< "Done face:" << Cf[facei]
                            << ", seeding outgoing edges" << endl;

                        // Seed all outgoing edges, i.e. those edges of
                        // which this face is owner
                        const labelList& fEdges = faceEdges[facei];
                        forAll(fEdges, fEdgei)
                        {
                            label myEdgei = fEdges[fEdgei];
                            if (edgeOwner[myEdgei] == facei)
                            {
                                if (isNewFrontEdge.set(edgei))
                                {
                                    Pout<< "Updating front edge:"
                                        << pp.localPoints()[e[0]]
                                        << pp.localPoints()[e[1]]
                                        << " with info from face:"
                                        << Cf[facei] << endl;

                                    newFrontEdges.append(myEdgei);
                                    allEdgeInfo[myEdgei] = allFaceInfo[facei];
                                }
                            }
                        }
                    }
                }
            }
        }


        if (returnReduce(newFrontEdges.empty(), sumOp<label>()))
        {
            break;
        }

        frontEdges.transfer(newFrontEdges);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
