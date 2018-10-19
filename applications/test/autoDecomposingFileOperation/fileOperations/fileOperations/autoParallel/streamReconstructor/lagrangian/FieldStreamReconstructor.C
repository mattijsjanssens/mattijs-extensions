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

\*---------------------------------------------------------------------------*/

#include "FieldStreamReconstructor.H"
#include "addToRunTimeSelectionTable.H"
#include "uFieldReconstructor.H"
#include "unallocatedIOPosition.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::FieldStreamReconstructor<Type>::reconstruct
(
    const IOobject& io,
    const bool,
    Ostream& os
) const
{
    const uFieldReconstructor& reconstructor =
        uFieldReconstructor::New(io.db());

    const PtrList<unallocatedFvMesh>& procMeshes = reconstructor.procMeshes();

    // Read field on proc meshes
    PtrList<cloud> procClouds(procMeshes.size());
    PtrList<IOField<Type>> procFields(procMeshes.size());

    label n = 0;
    forAll(procFields, proci)
    {
        Pout<< incrIndent;

        // Construct empty cloud
        procClouds.set
        (
            proci,
            new cloud
            (
                procMeshes[proci].thisDb(),
                "kinematicCloud"
            )
        );

        procFields.set
        (
            proci,
            new IOField<Type>
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    //io.local(),
                    procClouds[proci],
                    IOobject::MUST_READ,    //IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
        n += procFields[proci].size();

        Pout<< decrIndent;
    }

    IOField<Type> allField
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        n
    );

    n = 0;
    forAll(procFields, proci)
    {
        const IOField<Type>& procField = procFields[proci];
        SubList<Type>(allField, procField.size(), n) = procField;
        n += procField.size();
    }
    os << allField;
    return os.good();
}


template<class Type>
bool Foam::FieldStreamReconstructor<Type>::decompose
(
    const parUnallocatedFvFieldReconstructor& reconstructor,
    const unallocatedFvMesh& baseMesh,
    const IOobject& baseIO,

    const unallocatedFvMesh& thisMesh,
    const IOobject& thisIO,
    const bool,
    Ostream& os
) const
{
    Pout<< "*** FieldStreamReconstructor DEcomposing "
        << baseIO.objectPath() << endl;
    Pout<< "** FieldStreamReconstructor Decomposed "
        << baseIO.objectPath() << endl;
    return os.good();
}


template<class Type>
bool Foam::FieldStreamReconstructor<Type>::reconstruct
(
    const parUnallocatedFvFieldReconstructor& reconstructor,
    const regIOobject& thisIO,
    const bool,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    Pout<< "*** FieldStreamReconstructor reconstructin "
        << thisIO.objectPath() << endl;
    Pout<< "** FieldStreamReconstructor reconstruct "
        << thisIO.objectPath() << endl;
    return true;
}


// ************************************************************************* //
