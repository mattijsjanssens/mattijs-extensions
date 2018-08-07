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
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "zfp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

size_t compress(const scalarField& p, List<char>& buffer)
{
    zfp_type type = zfp_type_double;
    zfp_field* field = zfp_field_1d
    (
        reinterpret_cast<void*>(const_cast<scalarField&>(p).begin()),
        type,
        p.size()
    );

DebugVar(p.byteSize());


    /* allocate meta data for a compressed stream */
    zfp_stream* zfp = zfp_stream_open(nullptr);
    /* set compression mode and parameters  */
    zfp_stream_set_accuracy(zfp, 1e-6);
    //zfp_stream_set_precision(zfp, 32);

    /* allocate buffer for compressed data */
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
DebugVar(bufsize);
    buffer.setSize(bufsize);

    /* associate bit stream with allocated buffer */
    bitstream* stream = stream_open(buffer.begin(), bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    /* compress */
    size_t zfpsize = zfp_compress(zfp, field);
    buffer.setSize(zfpsize);
DebugVar(zfpsize);

    /* clean up */
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);

    return zfpsize;
}


void decompress(scalarField& p, List<char>& buffer)
{
    zfp_type type = zfp_type_double;
    zfp_field* field = zfp_field_1d
    (
        reinterpret_cast<void*>(p.begin()),
        type,
        p.size()
    );

    /* allocate meta data for a compressed stream */
    zfp_stream* zfp = zfp_stream_open(nullptr);
    /* set compression mode and parameters  */
    zfp_stream_set_accuracy(zfp, 1e-6);
    //zfp_stream_set_precision(zfp, 32);

    /* associate bit stream with allocated buffer */
    bitstream* stream = stream_open(buffer.begin(), buffer.byteSize());
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);
    if (!zfp_decompress(zfp, field))
    {
        FatalErrorInFunction << "Failed decompressing " << buffer.byteSize()
            << " buffer" << exit(FatalError);
    }

    /* clean up */
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    scalarField pData(10*p.size());
    label pDatai = 0;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        if (runTime.outputTime())
        {
DebugVar(pDatai);
DebugVar(pData.byteSize());
            {
                List<char> buf;
                compress(p.internalField(), buf);
                //DebugVar(buf);
                scalarField newP(p.internalField().size(), -123);
                decompress(newP, buf);
                forAll(newP, celli)
                {
                    Pout<< newP[celli] << '\t' << p[celli] << nl;
                }
            }


            if (pDatai == 10)
            {
                /* allocate meta data for the 2D array a[ny][nx] */
                zfp_type type = zfp_type_double;
                //size_t typesize = zfp_type_size(type);
                //zfp_field* field = zfp_field_alloc();
                //zfp_field_set_pointer(field, pData.begin());
                //zfp_field_set_type(field, zfp_type_double);
                //zfp_field_set_size_2d(field, 10, p.size());
                zfp_field* field = zfp_field_2d
                (
                    pData.begin(),
                    type,
                    p.size(),       // fastest varying
                    10
                );

                /* allocate meta data for a compressed stream */
                zfp_stream* zfp = zfp_stream_open(nullptr);
                /* set compression mode and parameters  */
                zfp_stream_set_accuracy(zfp, 1e-3);

                /* allocate buffer for compressed data */
                size_t bufsize = zfp_stream_maximum_size(zfp, field);
                void* buffer = malloc(bufsize);

                /* associate bit stream with allocated buffer */
                bitstream* stream = stream_open(buffer, bufsize);
                zfp_stream_set_bit_stream(zfp, stream);
                zfp_stream_rewind(zfp);

                /* compress */
                size_t zfpsize = zfp_compress(zfp, field);
DebugVar(zfpsize);
                fwrite(buffer, 1, zfpsize, stdout);

                /* clean up */
                zfp_field_free(field);
                zfp_stream_close(zfp);
                stream_close(stream);
                free(buffer);

                pDatai = 0;
            }
            else
            {
                // Shift up
                for (label i = 0; i < pDatai-1; i++)
                {
                    SubList<scalar>(pData, p.size(), pDatai*p.size()) =
                        SubList<scalar>(pData, p.size(), (pDatai+1)*p.size());
                }
                SubList<scalar>(pData, p.size(), pDatai*p.size()) =
                    p.internalField();

                pDatai++;
            }
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
