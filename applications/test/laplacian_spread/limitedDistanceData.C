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

#include "limitedDistanceData.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::limitedDistanceData<Type>& wDist
)
{
    return os
        << wDist.weights_ << token::SPACE
        << wDist.maxDist_ << token::SPACE
        << wDist.origin_ << token::SPACE
        << wDist.data_;
}


template<class Type>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::limitedDistanceData<Type>& wDist
)
{
    return is
        >> wDist.weights_ >> wDist.maxDist_
        >> wDist.origin_ >> wDist.data_;
}


// ************************************************************************* //
