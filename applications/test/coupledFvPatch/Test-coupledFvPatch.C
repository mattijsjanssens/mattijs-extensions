/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    test

Description
    CoupledFvPatch testing

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "argList.H"
#include "coupledFvPatch.H"
#include "coupledFvPatchFields.H"
#include "volFields.H"
#include "syncTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

const globalMeshData& gd = mesh.globalData();
{
    //const mapDistribute& edgeMap = gd.globalEdgeSlavesMap();

    labelField vals(identity(mesh.nEdges()));
    syncTools::syncEdgeList(mesh, vals, maxEqOp<label>(), labelMin);
    Pout<< "vals" << endl;
}

{
    const labelListList& pbf = gd.globalPointSlaves();
    const labelListList& pbtf = gd.globalPointTransformedSlaves();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const mapDistribute& map = gd.globalPointSlavesMap();

Pout<< "map.subMap():" << map.subMap() << endl;
Pout<< "map.constructMap():" << map.constructMap() << endl;


    pointField localPoints(cpp.localPoints());
    map.distribute(localPoints);


    forAll(pbf, pointI)
    {
        Pout<< "point:" << pointI << " at:" << localPoints[pointI]
            << endl;
        const labelList& pSlaves = pbf[pointI];
        forAll(pSlaves, i)
        {
            label cpPointI = pSlaves[i];
            Pout<< "cpPointI:" << cpPointI
                << " at:" << localPoints[cpPointI]
                << endl;
        }
        const labelList& pTrafoSlaves = pbtf[pointI];
        forAll(pTrafoSlaves, i)
        {
            label cpPointI = pTrafoSlaves[i];
            Pout<< "trafo cpPointI:" << cpPointI
                << " at:" << localPoints[cpPointI]
                << endl;
        }
        Pout<< endl;
    }
    return 0;
}

{
    const labelListList& pbf = gd.globalPointBoundaryFaces();
    const labelListList& pbtf = gd.globalPointTransformedBoundaryFaces();
    //const mapDistribute& map = gd.globalPointBoundaryFacesMap();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    forAll(pbf, pointI)
    {
        Pout<< "point:" << pointI << " at:" << cpp.localPoints()[pointI]
            << endl;
        const labelList& pFaces = pbf[pointI];
        forAll(pFaces, i)
        {
            label faceI = pFaces[i];
            Pout<< "face:" << faceI
                << " at:" << mesh.faceCentres()[faceI]
                << endl;
        }
        const labelList& pTrafoFaces = pbtf[pointI];
        forAll(pTrafoFaces, i)
        {
            label faceI = pTrafoFaces[i];
            Pout<< "trafo face:" << faceI
                << " at:" << mesh.faceCentres()[faceI]
                << endl;
        }
        Pout<< endl;
    }
}
return 0;

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const fvBoundaryMesh& fvbm = mesh.boundary();
    forAll(fvbm, patchI)
    {
        if (isA<coupledFvPatch>(fvbm[patchI]))
        {
            const coupledFvPatch& fvp =
                refCast<const coupledFvPatch>(fvbm[patchI]);

            Pout<< "Patch:" << fvp.name() << nl
                << incrIndent
                << indent << "Cf:" << fvp.Cf() << nl
                << indent << "Cn:" << fvp.Cn() << nl
                << indent << "Sf:" << fvp.Sf() << nl
                << indent << "magSf:" << fvp.magSf() << nl
                << indent << "nf:" << fvp.nf() << nl
                << indent << "delta:" << fvp.delta() << nl
                << indent << "weights:" << fvp.weights() << nl
                << indent << "deltaCoeffs:" << fvp.deltaCoeffs() << nl

                << indent << "forwardT:" << fvp.forwardT() << nl
                << indent << "reverseT:" << fvp.reverseT() << nl
                << decrIndent
                << endl;


            const coupledFvPatchScalarField& pf =
                refCast<const coupledFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                );

            Pout<< "PatchField:" << pf.type() << nl
                << incrIndent
                << indent << "value:" << pf << nl
                << indent << "snGrad:" << pf.snGrad() << nl
                << indent << "patchInternalField:"
                << pf.patchInternalField() << nl
                << indent << "patchNeighbourField:"
                << pf.patchNeighbourField() << nl
                << indent << "gradientInternalCoeffs:"
                << pf.gradientInternalCoeffs() << nl
                << indent << "gradientBoundaryCoeffs:"
                << pf.gradientBoundaryCoeffs() << nl
                << decrIndent
                << endl;
        }
    }

    Info<< "end" << endl;
}


// ************************************************************************* //
