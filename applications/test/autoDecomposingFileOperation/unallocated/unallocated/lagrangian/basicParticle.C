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

\*---------------------------------------------------------------------------*/

#include "basicParticle.H"
#include "IOstreams.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::basicParticle::sizeofPosition_
(
    offsetof(basicParticle, facei_) - offsetof(basicParticle, coordinates_)
);

const std::size_t Foam::basicParticle::sizeofFields_
(
    sizeof(basicParticle) - offsetof(basicParticle, coordinates_)
);

namespace Foam
{
    defineTypeNameAndDebug(basicParticle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicParticle::basicParticle
(
    Istream& is,
    bool readFields
)
:
    coordinates_(),
    celli_(-1),
    tetFacei_(-1),
    tetPti_(-1),
    facei_(-1),
    stepFraction_(0.0),
    origProc_(Pstream::myProcNo()),
    origId_(-1)
{
    if (is.format() == IOstream::ASCII)
    {
        is  >> coordinates_ >> celli_ >> tetFacei_ >> tetPti_;

        if (readFields)
        {
            is  >> facei_ >> stepFraction_ >> origProc_ >> origId_;
        }
    }
    else
    {
        if (readFields)
        {
            is.read(reinterpret_cast<char*>(&coordinates_), sizeofFields_);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&coordinates_), sizeofPosition_);
        }
    }

    // Check state of Istream
    is.check("basicParticle::basicParticle(Istream&, bool)");
}


Foam::basicParticle::basicParticle(const basicParticle& p, const label celli)
:
    coordinates_(p.coordinates_),
    celli_(celli),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::basicParticle::writePosition(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << coordinates_
            << token::SPACE << celli_
            << token::SPACE << tetFacei_
            << token::SPACE << tetPti_;
    }
    else
    {
        os.write(reinterpret_cast<const char*>(&coordinates_), sizeofPosition_);
    }

    // Check state of Ostream
    os.check("basicParticle::writePosition(Ostream& os, bool) const");
}


Foam::Ostream& Foam::operator<<(Ostream& os, const basicParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << p.coordinates_
            << token::SPACE << p.celli_
            << token::SPACE << p.tetFacei_
            << token::SPACE << p.tetPti_
            << token::SPACE << p.facei_
            << token::SPACE << p.stepFraction_
            << token::SPACE << p.origProc_
            << token::SPACE << p.origId_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&p.coordinates_),
            basicParticle::sizeofFields_
        );
    }

    return os;
}


// ************************************************************************* //
