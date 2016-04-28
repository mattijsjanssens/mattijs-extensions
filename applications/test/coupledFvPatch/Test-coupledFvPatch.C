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

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

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
                << indent
                << "Cf:" << fvp.Cf() << nl
                << "Cn:" << fvp.Cn() << nl
                << "Sf:" << fvp.Sf() << nl
                << "magSf:" << fvp.magSf() << nl
                << "nf:" << fvp.nf() << nl
                << "delta:" << fvp.delta() << nl
                << "weights:" << fvp.weights() << nl
                << "deltaCoeffs:" << fvp.deltaCoeffs() << nl

                << "forwardT:" << fvp.forwardT() << nl
                << "reverseT:" << fvp.reverseT() << nl
                << decrIndent
                << endl;


            const coupledFvPatchScalarField& pf =
                refCast<const coupledFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                );

            Pout<< "PatchField:" << pf.type() << nl
                << incrIndent
                << indent
                << "value:" << pf << nl
                << "snGrad:" << pf.snGrad() << nl
                << "patchInternalField:" << pf.patchInternalField() << nl
                << "patchNeighbourField:" << pf.patchNeighbourField() << nl
                << "gradientInternalCoeffs:" << pf.gradientInternalCoeffs()
                << nl
                << "gradientBoundaryCoeffs:" << pf.gradientBoundaryCoeffs()
                << nl
                << decrIndent
                << endl;
        }
    }

    Info<< "end" << endl;
}


// ************************************************************************* //
