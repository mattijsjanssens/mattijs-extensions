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
#include "SVD.H"

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


void Foam::fv::averageConstraint::average
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


void Foam::fv::averageConstraint::stencilWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
) const
{
    Pout<< "For sample:" << sample
        << " have ccs:" << donorCcs
        << endl;

    // Implicit least squares weighting
    // Number of donors
    label nD = donorCcs.size();

    weights.setSize(nD);

    // List for distance vectors and LSQ weights
    List<vector> d(nD);
    scalarList LSQw(nD);

    // Sum of weights
    scalar W = 0;

    // Sum of weighted distance vectors
    vector dw(Zero);

    RectangularMatrix<scalar> A(nD, 3);

    bool shortC = false;

    // Compute distance vectors and fill rectangular matrix
    forAll(donorCcs, j)
    {
        // Neighbour weights
        d[j] = donorCcs[j] - sample;

        // Check for short-circuiting if zero distance
        // is detected with respect to any donor
        if (mag(d[j]) < ROOTVSMALL)
        {
            shortC = true;
            break;
        }

        LSQw[j] = 1.0/magSqr(d[j]);

        // T matrix
        vector wd = LSQw[j]*d[j];
        A[j][0] = wd.x();
        A[j][1] = wd.y();
        A[j][2] = wd.z();

        // Sum of weighted distance vectors
        dw += wd;

        // Sum of weights
        W += LSQw[j];
    }

    if (!shortC)
    {
        // Use Singular Value Decomposition to avoid problems
        // with 1D, 2D stencils
        SVD svd(A.T()*A, SMALL);

        // Least squares vectors
        RectangularMatrix<scalar> ATAinvAT(svd.VSinvUt()*A.T());

        scalar saveDiag(W);

        // Diagonal coefficient
        for (label i = 0; i < 3; i++)
        {
            // Get row
            scalarList Row(UList<scalar>(ATAinvAT[i], nD));

            forAll(donorCcs, j)
            {
                W -= Row[j]*LSQw[j]*dw[i];
            }
        }

        if (mag(W) < SMALL*mag(saveDiag))
        {
            shortC = true;
        }
        else
        {
            // Compute final neighbour weights with  additional scaling
            forAll(donorCcs, j)
            {
                weights[j] =
                (
                    LSQw[j]
                - ATAinvAT[0][j]*LSQw[j]*dw[0]
                - ATAinvAT[1][j]*LSQw[j]*dw[1]
                - ATAinvAT[2][j]*LSQw[j]*dw[2]
                )
            /W;
            }
        }
    }

    if (shortC)
    {
        // Matrix ill conditioned. Use straight injection from central
        // donor.
        weights = 0.0;
        weights[0] = 1.0;
    }
}
void Foam::fv::averageConstraint::interpolate
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
    const auto& C = mesh_.C();


    // Collect neighbour cell centres

    DynamicList<point> donorCcs;
    DynamicList<label> faceIDs;
    DynamicList<label> patchIDs;
    for (const label celli : cells_)
    {
        Pout<< "OLD: For cell:" << celli
            << " at:" << mesh_.cellCentres()[celli]
            << " diag:" << eqn.diag()[celli]
            << " source:" << eqn.source()[celli]
            << endl;

        // Count number of valid contributions (most boundary faces will not
        // add a contribution - only coupled ones will)
        donorCcs.clear();
        faceIDs.clear();
        patchIDs.clear();
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
                donorCcs.append(C[nbrCelli]);
                faceIDs.append(facei);
                patchIDs.append(-1);
            }
            else
            {
                const label patchi = mesh_.boundaryMesh().whichPatch(facei);
                const auto& pp = mesh_.boundaryMesh()[patchi];

                if (eqn.boundaryCoeffs()[patchi].size())
                {
                    const label patchFacei = pp.whichFace(facei);
                    donorCcs.append(C.boundaryField()[patchi][patchFacei]);
                    faceIDs.append(facei);
                    patchIDs.append(patchi);
                }
            }
        }
        scalarList weights;
        stencilWeights
        (
            C[celli],
            donorCcs,
            weights
        );

        const scalar sumWeights = sum(weights);

//XXXXXX
        // Insert as averaging stencil
        const scalar oldDiag = eqn.diag()[celli];

        eqn.source()[celli] = 0.0;
        eqn.diag()[celli]*= sumWeights;

        // And set psi already to average
        eqn.psi()[celli] = Zero;

        // Set all off diagonal contributions
        forAll(faceIDs, i)
        {
            const scalar w = oldDiag*weights[i];
            const label facei = faceIDs[i];

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
                    eqn.upper()[facei] = w;
                }
                else
                {
                    if (celli == own[facei])
                    {
                        eqn.upper()[facei] = w;
                    }
                    else
                    {
                        eqn.lower()[facei] = w;
                    }
                }

                eqn.psi()[celli] += w*eqn.psi()[nbrCelli];
            }
            else
            {
                const label patchi = mesh_.boundaryMesh().whichPatch(facei);
                const auto& pp = mesh_.boundaryMesh()[patchi];

                if (eqn.boundaryCoeffs()[patchi].size())
                {
                    const label patchFacei = pp.whichFace(facei);
                    eqn.boundaryCoeffs()[patchi][patchFacei] = w;
                    eqn.psi()[celli] +=
                        w*eqn.psi().boundaryField()[patchi][patchFacei];
                }
            }
        }
        eqn.psi()[celli] /= sumWeights;

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


void Foam::fv::averageConstraint::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    //average(eqn, field);
    interpolate(eqn, fieldi);
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
