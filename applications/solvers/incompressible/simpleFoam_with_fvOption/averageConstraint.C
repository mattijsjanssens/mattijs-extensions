/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 M. Janssens
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

#include "averageConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(averageConstraint, 0);
        addToRunTimeSelectionTable
        (
            option,
            averageConstraint,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::averageConstraint::averageConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::averageConstraint::print
(
    fvMatrix<scalar>& eqn,
    const label celli
) const
{
    const cellList& cells = mesh_.cells();
    const labelUList& own = mesh_.owner();
    //const labelUList& nei = mesh_.neighbour();

    Pout<< "For cell:" << celli
        << " at:" << mesh_.cellCentres()[celli]
        << " value:" << eqn.psi()[celli]
        << " diag:" << eqn.diag()[celli]
        << " source:" << eqn.source()[celli]
        << endl;


    for (const label facei : cells[celli])
    {
        if (mesh_.isInternalFace(facei))
        {
            //const label nbrCelli
            //(
            //    (celli == own[facei])
            //  ? nei[facei]
            //  : own[facei]
            //);

            if (eqn.symmetric())
            {
                // Only upper is set
                Pout<< "    internal face:" << facei
                    << " at:" << mesh_.faceCentres()[facei]
                    << " upper:" << eqn.upper()[facei] << endl;
            }
            else
            {
                if (celli == own[facei])
                {
                    Pout<< "    internal face:" << facei
                        << " at:" << mesh_.faceCentres()[facei]
                        << " upper:" << eqn.upper()[facei] << endl;
                }
                else
                {
                    Pout<< "    internal face:" << facei
                        << " at:" << mesh_.faceCentres()[facei]
                        << " lower:" << eqn.lower()[facei] << endl;
                }
            }
        }
        else
        {
            const label patchi = mesh_.boundaryMesh().whichPatch(facei);

            if (eqn.boundaryCoeffs()[patchi].size())
            {
                const label patchFacei =
                    mesh_.boundaryMesh()[patchi].whichFace(facei);

                Pout<< "    boundary face:" << facei
                    << " at:" << mesh_.faceCentres()[facei]
                    << " boundaryCoeffs:"
                    << eqn.boundaryCoeffs()[patchi][patchFacei]
                    << " internalCoeffs:"
                    << eqn.internalCoeffs()[patchi][patchFacei]
                    << endl;

                // What to do about neighbour value?
            }
        }
    }
}


void Foam::fv::averageConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    Pout
        << "averageConstraint::constrain"
        << " for field " << eqn.psi().name() << endl;

    if (!eqn.symmetric() && !eqn.asymmetric())
    {
        return;
    }

    const cellList& cells = mesh_.cells();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();


    // Keep diagonal whatever it is and set off-diagonals

    for (const label celli : cells_)
    {
        Pout<< "OLD: For cell:" << celli
            << " at:" << mesh_.cellCentres()[celli]
            << " diag:" << eqn.diag()[celli]
            << " source:" << eqn.source()[celli]
            << endl;

        // Count number of valid contributions (most boundary faces will not
        // add a contribution - only coupled ones will)
        label nFaces = 0;
        for (const label facei : cells[celli])
        {
            if (mesh_.isInternalFace(facei))
            {
                nFaces++;
            }
            else
            {
                const label patchi = mesh_.boundaryMesh().whichPatch(facei);

                if (eqn.boundaryCoeffs()[patchi].size())
                {
                    nFaces++;
                }
            }
        }

        eqn.source()[celli] = 0.0;


        const scalar offDiag = -eqn.diag()[celli]/nFaces;
        eqn.psi()[celli] = Zero;


        // Set all off diagonal contributions
        for (const label facei : cells[celli])
        {
            if (mesh_.isInternalFace(facei))
            {
                const label nbrCelli
                (
                    (celli == own[facei])
                  ? nei[facei]
                  : own[facei]
                );

                if (eqn.symmetric())
                {
                    // Only upper is set
                    eqn.upper()[facei] = offDiag;
                }
                else
                {
                    if (celli == own[facei])
                    {
                        eqn.upper()[facei] = offDiag;
                    }
                    else
                    {
                        eqn.lower()[facei] = offDiag;
                    }
                }

                eqn.psi()[celli] += eqn.psi()[nbrCelli]/nFaces;
            }
            else
            {
                const label patchi = mesh_.boundaryMesh().whichPatch(facei);

                if (eqn.boundaryCoeffs()[patchi].size())
                {
                    const label patchFacei =
                        mesh_.boundaryMesh()[patchi].whichFace(facei);
                    eqn.boundaryCoeffs()[patchi][patchFacei] = offDiag;

                    // What to do about neighbour value?
                }
            }
        }

        {
            print(eqn, celli);
            const labelList& cCells = mesh_.cellCells()[celli];
            for (const label celli : cCells)
            {
                Pout<< "    Neighbour:" << celli
                    << " at:" << mesh_.cellCentres()[celli]
                    << " value:" << eqn.psi()[celli]
                    << endl;
            }
        }
    }
}


bool Foam::fv::averageConstraint::read(const dictionary& dict)
{
    Pout<< "averageConstraint::read dict:" << dict << endl;

    if (fv::cellSetOption::read(dict))
    {
        fieldNames_ = dict.get<wordList>("fields");
DebugVar(fieldNames_);
        fv::option::resetApplied();

        return true;
    }

    return false;
}


// ************************************************************************* //
