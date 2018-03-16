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
#include "PatchEdgeFaceWave.H"
#include "patchIntegrateInfo.H"
#include "EdgeMap.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class PatchType>
tmp<Field<T>> fecLinearInterpolate
(
    const polyMesh& mesh,
    const PatchType& pp,
    const labelUList& meshPoints,
    const scalarField& faceCentres,
    const scalarField& edgeCentres,
    const Field<T>& faceFld
)
{
    if (edgeCentres.size() != pp.nEdges() || faceFld.size() != pp.size())
    {
        FatalErrorInFunction << "Incorrect sizes" << exit(FatalError);
    }

    tmp<Field<T>> tsumFld(new Field<T>(pp.nEdges(), 0));
    Field<T>& sumFld = tsumFld.ref();
    scalarField sumWeights(pp.nEdges(), 0);

    // Sum local face contributions
    {
        const labelListList& edgeFaces = pp.edgeFaces();
        forAll(edgeFaces, edgei)
        {
            const labelList& eFaces = edgeFaces[edgei];

            forAll(eFaces, eFacei)
            {
                label facei = eFaces[eFacei];

                scalar d = mag(faceCentres[facei]-edgeCentres[edgei]);
                scalar w = 1.0/max(SMALL, d);
                sumFld[edgei] += w*faceFld[facei];
                sumWeights[edgei] += w;
            }
        }
    }


    // Sum coupled contributions
    {
        const globalMeshData& gd = mesh.globalData();
        const indirectPrimitivePatch& cpp = gd.coupledPatch();

        // Build map from pp edges (in mesh indices) to edge indices
        EdgeMap<label> edgeToIndex(2*pp.nEdges());
        {
            const edgeList& edges = pp.edges();
            forAll(edges, edgei)
            {
                const edge& e = edges[edgei];
                const edge meshE(meshPoints[e[0]], meshPoints[e[1]]);
                edgeToIndex.insert(meshE, edgei);
            }
        }

        // Fill coupled edge data
        Field<T> cppSum(cpp.nEdges(), 0.0);
        scalarField cppWeight(cpp.nEdges(), 0.0);
        {
            const edgeList& cppEdges = cpp.edges();
            const labelList& mp = cpp.meshPoints();
            forAll(cppEdges, cppEdgei)
            {
                const edge& e = cppEdges[cppEdgei];
                const edge meshE(mp[e[0]], mp[e[1]]);
                EdgeMap<label>::const_iterator iter = edgeToIndex.find(meshE);
                if (iter != edgeToIndex.end())
                {
                    cppSum[cppEdgei] = sumFld[iter()];
                    cppWeight[cppEdgei] = sumWeights[iter()];
                }
            }
        }

        // Add contributions
        globalMeshData::syncData
        (
            cppSum,
            gd.globalEdgeSlaves(),
            gd.globalEdgeTransformedSlaves(),
            gd.globalEdgeSlavesMap(),
            gd.globalTransforms(),
            plusEqOp<T>(),
            mapDistribute::transform()
        );
        globalMeshData::syncData
        (
            cppWeight,
            gd.globalEdgeSlaves(),
            gd.globalEdgeTransformedSlaves(),
            gd.globalEdgeSlavesMap(),
            gd.globalTransforms(),
            plusEqOp<scalar>(),
            mapDistribute::transform()
        );

        // Extract
        {
            const edgeList& cppEdges = cpp.edges();
            const labelList& mp = cpp.meshPoints();
            forAll(cppEdges, cppEdgei)
            {
                const edge& e = cppEdges[cppEdgei];
                const edge meshE(mp[e[0]], mp[e[1]]);
                EdgeMap<label>::const_iterator iter = edgeToIndex.find(meshE);
                if (iter != edgeToIndex.end())
                {
                    sumFld[iter()] = cppSum[cppEdgei];
                    sumWeights[iter()] = cppWeight[cppEdgei];
                }
            }
        }
    }


    // Normalise
    forAll(sumWeights, edgei)
    {
        if (sumWeights[edgei] > VSMALL)
        {
            sumFld[edgei] /= sumWeights[edgei];
        }
    }

    return tsumFld;
}



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
    //const vectorField& Cf = fvp.Cf();
    const polyPatch& pp = fvp.patch();
    //const labelListList& edgeFaces = pp.edgeFaces();
    //const labelListList& faceEdges = pp.faceEdges();

//// Addressing: owner = geometrically lower face, neighbour = geometrically
//// higher face
//labelList edgeOwner(pp.nEdges(), -1);
//labelList edgeNeighbour(pp.nEdges(), -1);
//labelList nInEdges(pp.size(), 0);
//{
//    forAll(edgeFaces, edgei)
//    {
//        const edge& e = pp.edges()[edgei];
//        const point& p0 = pp.localPoints()[e[0]];
//        const point& p1 = pp.localPoints()[e[1]];
//        point edgeMid(0.5*(p0+p1));
//
//        const labelList& eFaces = edgeFaces[edgei];
//        forAll(eFaces, i)
//        {
//            label facei = eFaces[i];
//            const point& fc = Cf[facei];
//
//            if (fc.x() < edgeMid.x())
//            {
//                edgeOwner[edgei] = facei;
//            }
//            else
//            {
//                edgeNeighbour[edgei] = facei;
//
//                if (eFaces.size() > 1)
//                {
//                    nInEdges[facei]++;
//                }
//            }
//        }
//    }
//}


    // Get component to integrate along
    const scalarField faceCentres(pp.faceCentres().component(0));
    scalarField edgeCentres(pp.nEdges());
    {
        forAll(pp.edges(), edgei)
        {
            const edge& e = pp.edges()[edgei];
            const point& p0 = mesh.points()[pp.meshPoints()[e[0]]];
            const point& p1 = mesh.points()[pp.meshPoints()[e[1]]];
            edgeCentres[edgei] = 0.5*(p0[0]+p1[0]);
        }
    }


// // Test interpolation
// {
//     tmp<scalarField> tedgeFld
//     (
//         fecLinearInterpolate
//         (
//             pp,
//             faceCentres,
//             edgeCentres,
//             faceCentres
//         )
//     );
//
//     forAll(pp.edges(), edgei)
//     {
//         //const edge& e = pp.edges()[edgei];
//         //const point& p0 = mesh.points()[pp.meshPoints()[e[0]]];
//         //const point& p1 = mesh.points()[pp.meshPoints()[e[1]]];
//
//         const labelList& eFaces = pp.edgeFaces()[edgei];
//
//         Pout<< "At edge:" << edgeCentres[edgei]
//             << " have face values:"
//             << UIndirectList<scalar>(faceCentres, eFaces)
//             << " and edge value:" << tedgeFld()[edgei] << endl;
//     }
//     return 0;
// }





    List<patchIntegrateInfo> allEdgeInfo(pp.nEdges());
    List<patchIntegrateInfo> allFaceInfo(pp.size());

    DynamicList<label> changedEdges;
    DynamicList<patchIntegrateInfo> changedInfo;
    {
        // Seed
        forAll(pp.edges(), edgei)
        {
            if (edgeCentres[edgei] < 1e-5)
            {
                const edge& e = pp.edges()[edgei];
                Pout<< "Front edge:"
                    << pp.localPoints()[e[0]]
                    << pp.localPoints()[e[1]]
                    << endl;

                changedEdges.append(edgei);
                changedInfo.append(patchIntegrateInfo(0.0));
            }
        }
    }


    scalarField leftFaceMask(pos(0.05-pp.faceCentres().component(0)));
    DebugVar(leftFaceMask);

    // scalarField leftEdgeMask(pp.nEdges());
    // {
    //     forAll(pp.edges(), edgei)
    //     {
    //         const edge& e = pp.edges()[edgei];
    //         const point& p0 = mesh.points()[pp.meshPoints()[e[0]]];
    //         const point& p1 = mesh.points()[pp.meshPoints()[e[1]]];
    //         const point eMid(0.5*(p0+p1));
    //         if (eMid[0] < 0.05)
    //         {
    //             leftEdgeMask[edgei] = 1.0;
    //         }
    //         else
    //         {
    //             leftEdgeMask[edgei] = 0.0;
    //         }
    //     }
    // }
    // DebugVar(leftEdgeMask);


    const scalarField faceDensity(leftFaceMask*scalarField(pp.size(), 1.0));
    const scalarField edgeDensity
    (
        fecLinearInterpolate
        (
            mesh,
            pp,
            pp.meshPoints(),
            faceCentres,
            edgeCentres,
            faceDensity
        )
    );

    patchIntegrateInfo::trackData td
    (
        faceCentres,
        edgeCentres,
        faceDensity,
        edgeDensity
    );


    // Walk
    PatchEdgeFaceWave
    <
        primitivePatch,
        patchIntegrateInfo,
        patchIntegrateInfo::trackData
    > calc
    (
        mesh,
        pp,
        changedEdges,
        changedInfo,
        allEdgeInfo,
        allFaceInfo,
        returnReduce(pp.nEdges(), sumOp<label>()),
        td
    );


DebugVar(allEdgeInfo);
DebugVar(allFaceInfo);


    // Extract as patchField
    volScalarField vsf
    (
        IOobject
        (
            "patchDist",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("patchDist", dimLength, 0.0)
    );
    scalarField pf(vsf.boundaryField()[pp.index()].size());
    forAll(pf, facei)
    {
        pf[facei] = allFaceInfo[facei].value();
    }
    vsf.boundaryFieldRef()[pp.index()] = pf;

    Info<< "Writing patchDist volScalarField to " << runTime.value()
        << endl;

    vsf.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
