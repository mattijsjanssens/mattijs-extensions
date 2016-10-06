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

Description
    Wall plane patch

\*---------------------------------------------------------------------------*/

#include "fixedConstraintPointPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(fixedConstraintPointPatch, 0);

addToRunTimeSelectionTable
(
    pointPatch,
    fixedConstraintPointPatch,
    dictionary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedConstraintPointPatch::fixedConstraintPointPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const pointBoundaryMesh& bm,
    const word& patchType
)
:
    pointPatch(bm),
    name_(name),
    index_(index),
    meshPoints_(dict.lookup("meshPoints")),
    constraints_(dict.lookup("constraints"))
{
    if (meshPoints_.size() != constraints_.size())
    {
        FatalErrorInFunction << "patch " << name_
            << " size of meshPoints " << meshPoints_.size()
            << " differs from size of constraints " << constraints_.size()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::fixedConstraintPointPatch::localPoints() const
{
    if (!localPointsPtr_.valid())
    {
        const pointField& points = boundaryMesh().mesh()().points();

        localPointsPtr_.reset(new pointField(points, meshPoints()));
    }
    return localPointsPtr_();
}


const Foam::vectorField& Foam::fixedConstraintPointPatch::pointNormals() const
{
    if (!pointNormalsPtr_.valid())
    {
        pointNormalsPtr_.reset(new vectorField(constraints_.size()));
        vectorField& pointNormals = pointNormalsPtr_();
        forAll(constraints_, i)
        {
            pointNormals[i] = constraints_[i].second();
        }
    }
    return pointNormalsPtr_();
}


void Foam::fixedConstraintPointPatch::setConstraints
(
    const List<pointConstraint>& pc
)
{
    constraints_ = pc;
}


void Foam::fixedConstraintPointPatch::write(Ostream& os) const
{
    pointPatch::write(os);
    os.writeKeyword("meshPoints")
        << meshPoints() << token::END_STATEMENT << nl;
    os.writeKeyword("constraints")
        << constraints() << token::END_STATEMENT << nl;
}


// ************************************************************************* //
