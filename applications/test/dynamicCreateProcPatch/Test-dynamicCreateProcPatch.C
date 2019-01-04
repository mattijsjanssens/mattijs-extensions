/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Test-fieldMapping

Description
    Test app for mapping of fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
//#include "zeroGradientFvPatchFields.H"
#include "processorFvPatchFields.H"
#include "fvMeshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// bool notEqual(const scalar s1, const scalar s2, const scalar tol)
// {
//     return mag(s1-s2) > tol;
// }

template<class GeoField>
void savePatchValues
(
    objectRegistry& values,
    const objectRegistry& db,
    const label patchi
)
{
    typedef IOField<typename GeoField::value_type> ValType;

    const wordList fldNames(db.sortedNames(GeoField::typeName));

    forAll(fldNames, i)
    {
        const GeoField& fld = db.lookupObject<GeoField>(fldNames[i]);
        const typename GeoField::Patch& pfld = fld.boundaryField()[patchi];

        IOobject io(fld, values);
        io.readOpt() = IOobject::NO_READ;
        io.writeOpt() = IOobject::NO_WRITE;
        io.registerObject() = true;

        regIOobject::store(new ValType(io, pfld));
    }
}


template<class GeoField>
void setPatchValues
(
    objectRegistry& db,
    const objectRegistry& values,
    const mapPolyMesh& mpm,
    const label fromPatchi,
    const label patchi
)
{
    typedef IOField<typename GeoField::value_type> ValType;

    const labelList& oldToNew = mpm.reverseFaceMap();
    const label oldStart = mpm.oldPatchStarts()[fromPatchi];
    const label oldSize = mpm.oldPatchSizes()[fromPatchi];
    const label newStart = mpm.mesh().boundaryMesh()[patchi].start();


    const wordList fldNames(db.sortedNames(GeoField::typeName));

    forAll(fldNames, i)
    {
        GeoField& fld = db.lookupObjectRef<GeoField>(fldNames[i]);
        typename GeoField::Patch& pfld = fld.boundaryFieldRef()[patchi];

        // Look up values
        const ValType& vals = values.lookupObject<ValType>(fldNames[i]);

        // Start off with copy
        Field<typename GeoField::value_type> fvp(pfld);
        for
        (
            label oldFacei = oldStart;
            oldFacei < oldStart+oldSize;
            oldFacei++
        )
        {
            label oldPatchFacei = oldFacei-oldStart;
            label procFacei = oldToNew[oldFacei]-newStart;
            fvp[procFacei] = vals[oldPatchFacei];
        }

        pfld = fvp;
    }
}


// Main program:

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const faceZoneMesh& fzm = mesh.faceZones();

    const label patchi = pbm.findPatchID("matchPatch");
    if (patchi == -1)
    {
        FatalErrorInFunction << "no patch matchPatch" << exit(FatalError);
    }

    if (Pstream::nProcs() != 2)
    {
        FatalErrorInFunction << "run on two processors only"
            << exit(FatalError);
    }


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


    // Add processor patch between 0 and 1
    label procPatchi = -1;
    {
        label nbrProci = 1-Pstream::myProcNo();

        processorPolyPatch pp
        (
            0,              // size
            mesh.nFaces(),
            pbm.size(),
            pbm,
            Pstream::myProcNo(),    // my processor
            nbrProci                // nbr processor
        );

        procPatchi = fvMeshTools::addPatch
        (
            mesh,
            pp,
            dictionary(),   // optional per field patchField
            processorFvPatchField<scalar>::typeName,
            false           // not parallel sync
        );
    }



    {
        polyTopoChange meshMod(mesh);

        // Save old patch values
        //const scalarField vals(ccX.boundaryField()[patchi]);
        objectRegistry obr(IOobject("boundaryValues", mesh));

        savePatchValues<volScalarField>(obr, mesh, patchi);
        savePatchValues<volVectorField>(obr, mesh, patchi);


        // Move faces from patch
        const polyPatch& pp = pbm[patchi];
        forAll(pp, i)
        {
            label facei = pp.start()+i;
            
            const label zonei = fzm.whichZone(facei);
            bool flip = false;
            if (zonei != -1)
            {
                flip = fzm[zonei].flipMap()[fzm[zonei].whichFace(facei)];
            }

            meshMod.modifyFace
            (
                mesh.faces()[facei],
                facei,
                mesh.faceOwner()[facei],
                -1,
                false,          // no need to flip face
                procPatchi,
                zonei,
                flip
            );
        }

        // Change mesh and inflate
        Info<< "Actually changing mesh" << nl << endl;
        autoPtr<mapPolyMesh> mpm = meshMod.changeMesh(mesh, false);

        Info<< "Mapping fields" << nl << endl;
        mesh.updateMesh(mpm);

        setPatchValues<volScalarField>(mesh, obr, mpm(), patchi, procPatchi);
        setPatchValues<volVectorField>(mesh, obr, mpm(), patchi, procPatchi);


        // Move mesh (since morphing does not do this)
        if (mpm().hasMotionPoints())
        {
            Info<< "Moving mesh" << nl << endl;
            mesh.movePoints(mpm().preMotionPoints());
        }
    }

    runTime++;

    Info<< "Writing mesh and fields" << nl << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
