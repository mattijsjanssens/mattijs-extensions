/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 M. Janssens
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

#include "blockConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(blockConstraint, 0);
        addToRunTimeSelectionTable
        (
            option,
            blockConstraint,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::blockConstraint::blockConstraint
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

//void Foam::fv::blockConstraint::print
//(
//    fvMatrix<scalar>& eqn,
//    const label celli
//) const
//{
//    const cellList& cells = mesh_.cells();
//    const labelUList& own = mesh_.owner();
//    //const labelUList& nei = mesh_.neighbour();
//
//    Pout<< "For cell:" << celli
//        << " at:" << mesh_.cellCentres()[celli]
//        << " value:" << eqn.psi()[celli]
//        << " diag:" << eqn.diag()[celli]
//        << " source:" << eqn.source()[celli]
//        << endl;
//
//
//    for (const label facei : cells[celli])
//    {
//        if (mesh_.isInternalFace(facei))
//        {
//            //const label nbrCelli
//            //(
//            //    (celli == own[facei])
//            //  ? nei[facei]
//            //  : own[facei]
//            //);
//
//            if (eqn.symmetric())
//            {
//                // Only upper is set
//                Pout<< "    internal face:" << facei
//                    << " at:" << mesh_.faceCentres()[facei]
//                    << " upper:" << eqn.upper()[facei] << endl;
//            }
//            else
//            {
//                if (celli == own[facei])
//                {
//                    Pout<< "    internal face:" << facei
//                        << " at:" << mesh_.faceCentres()[facei]
//                        << " upper:" << eqn.upper()[facei] << endl;
//                }
//                else
//                {
//                    Pout<< "    internal face:" << facei
//                        << " at:" << mesh_.faceCentres()[facei]
//                        << " lower:" << eqn.lower()[facei] << endl;
//                }
//            }
//        }
//        else
//        {
//            const label patchi = mesh_.boundaryMesh().whichPatch(facei);
//
//            if (eqn.boundaryCoeffs()[patchi].size())
//            {
//                const label patchFacei =
//                    mesh_.boundaryMesh()[patchi].whichFace(facei);
//
//                Pout<< "    boundary face:" << facei
//                    << " at:" << mesh_.faceCentres()[facei]
//                    << " boundaryCoeffs:"
//                    << eqn.boundaryCoeffs()[patchi][patchFacei]
//                    << " internalCoeffs:"
//                    << eqn.internalCoeffs()[patchi][patchFacei]
//                    << endl;
//
//                // What to do about neighbour value?
//            }
//        }
//    }
//}


void Foam::fv::blockConstraint::interpolate
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    Pout
        << "blockConstraint::constrain"
        << " for field " << eqn.psi().name() << endl;

    if (!eqn.symmetric() && !eqn.asymmetric())
    {
        return;
    }

    const cellList& cells = mesh_.cells();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    bitSet isSetCell(mesh_.nCells());
    isSetCell.set(cells_);

    bitSet isSetFace(mesh_.nBoundaryFaces());
    for (label bFacei = 0; bFacei < mesh_.nBoundaryFaces(); bFacei++)
    {
        const label own = mesh_.faceOwner()[mesh_.nInternalFaces()+bFacei];
        if (isSetCell[own])
        {
            isSetFace.set(bFacei);
        }
    }
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        isSetFace,
        orEqOp<unsigned int>()
    );


    // Set all coefficients to external cells to zero

    for (const label celli : cells_)
    {
        Pout<< "OLD: For cell:" << celli
            << " at:" << mesh_.cellCentres()[celli]
            << " diag:" << eqn.diag()[celli]
            << " source:" << eqn.source()[celli]
            << endl;

        for (const label facei : cells[celli])
        {
            if (mesh_.isInternalFace(facei))
            {
                if (isSetCell[own[facei]] && !isSetCell[nei[facei]])
                {
                    // Cut connection to neighbour
                    eqn.lower()[facei] = 0.0;
                }
                else if (!isSetCell[own[facei]] && isSetCell[nei[facei]])
                {
                    // Cut connection to neighbour
                    eqn.upper()[facei] = 0.0;
                }
            }
            else
            {
                const label patchi = mesh_.boundaryMesh().whichPatch(facei);
                const auto& pp = mesh_.boundaryMesh()[patchi];

                if (eqn.boundaryCoeffs()[patchi].size())
                {
                    const label patchFacei = pp.whichFace(facei);
                    const label bFacei = facei-mesh_.nInternalFaces();

                    if (!isSetCell[own[facei]] && isSetFace[bFacei])
                    {
                        // Cut connection to neighbour
                        eqn.boundaryCoeffs()[patchi][patchFacei] = 0.0;
                    }
                }
            }
        }
    }
}


void Foam::fv::blockConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    interpolate(eqn, fieldi);
}


bool Foam::fv::blockConstraint::read(const dictionary& dict)
{
    Pout<< "blockConstraint::read dict:" << dict << endl;

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
