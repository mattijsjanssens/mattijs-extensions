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
#include "correctedSnGrad.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

//     Info<< "Reading field p\n" << endl;
//     volScalarField p
//     (
//         IOobject
//         (
//             "p",
//             runTime.timeName(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh
//     );
    Info<< "Constructing field p\n" << endl;
    volScalarField p
    (
        "p",
        mesh.C().component(vector::X)
    );

    const surfaceScalarField& deltaCoeffs =
        fv::correctedSnGrad<scalar>(mesh).deltaCoeffs(p);

    const fvBoundaryMesh& fvbm = mesh.boundary();
    forAll(fvbm, patchI)
    {
        if (isA<coupledFvPatch>(fvbm[patchI]))
        {
            const coupledFvPatch& fvp =
                refCast<const coupledFvPatch>(fvbm[patchI]);

            Pout<< "Patch:" << fvp.name() << " type:" << fvp.type() << nl
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


//             const coupledFvPatchScalarField& pf =
//                 refCast<const coupledFvPatchScalarField>
//                 (
//                     p.boundaryField()[patchI]
//                 );
            const fvPatchScalarField& pf = p.boundaryField()[patchI];

            Pout<< "PatchField:" << pf.type() << endl;

            const fvsPatchScalarField& pDeltaCoeffs =
            deltaCoeffs.boundaryField()[patchI];


            Pout<< incrIndent
                << indent << "value:" << static_cast<const scalarField&>(pf)
                << endl;

             Pout<< indent << "snGrad:" << pf.snGrad(pDeltaCoeffs) << nl
                 << indent << "patchInternalField:"
                 << pf.patchInternalField() << nl
                 << indent << "patchNeighbourField:"
                 << pf.patchNeighbourField() << endl;
 
 
             Pout<< indent << "gradientInternalCoeffs:"
                 << pf.gradientInternalCoeffs(pDeltaCoeffs) << nl
                 << indent << "gradientBoundaryCoeffs:"
                 << pf.gradientBoundaryCoeffs(pDeltaCoeffs) << nl
                 << decrIndent
                 << endl;
        }
    }

    volVectorField gradP("gradP", fvc::grad(p));
    Pout<< "gradP:" << gradP << endl;
    Info<< "Writing " << gradP.name() << " to " << runTime.timeName() << endl;
    gradP.write();


    runTime++;
    Info<< "Writing p to " << runTime.timeName() << endl;
    p.correctBoundaryConditions();
    p.write();


    Info<< "end" << endl;
}


// ************************************************************************* //
