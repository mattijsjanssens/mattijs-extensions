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

//#include "fvFieldDecomposer.H"
#include "parUnallocatedFvFieldReconstructor.H"
#include "uVolFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoField>
bool Foam::fileOperations::autoParallelFileOperation::reconstructAndWrite
(
    const GeoField& procFld,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    Pout<< "**REconstructing " << procFld.name() << endl;

    if (!reconstructorPtr_.valid())
    {
        const unallocatedFvMesh& baseM = baseMesh(procFld.time());

        Info<< "Creating reconstructor" << nl << endl;
        reconstructorPtr_ = new parUnallocatedFvFieldReconstructor
        (
            baseM,
            procFld.mesh(),
            distMapPtr_()
        );
    }
    const parUnallocatedFvFieldReconstructor& reconstructor =
        reconstructorPtr_();

    // Map local field onto baseMesh
    tmp<uVolScalarField> tfld
    (
        reconstructor.reconstructFvVolumeField(procFld)
    );

    // Write master field to parent
    bool state = true;
    if (Pstream::master())
    {
        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;
        {
            Pout<< "**Writign " << tfld().objectPath() << endl;
            //state = tfld().write();
            state = tfld().writeObject(fmt, ver, cmp, true);
        }
        Pstream::parRun() = oldParRun;
    }

    return state;
}


// ************************************************************************* //
