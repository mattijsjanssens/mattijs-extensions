/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "parProfiling.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamGlobals.H"
#include "PstreamReduceOps.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(parProfiling, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        parProfiling,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::parProfiling::parProfiling
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::parProfiling::~parProfiling()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::parProfiling::execute()
{
    scalarList times(Pstream::nProcs());
    times[Pstream::myProcNo()] = PstreamGlobals::reduceTime_;
    Pstream::gatherList(times);

//     Info<< "Reduce : min/max/avg : " << reduceMin << '/' << reduceMax
//         << '/' << reduceSum/Pstream::nProcs() << nl;

    DebugVar(times);

//
// std::cout<< "** Timings" << std::endl
//     << "\treduce   : " << PstreamGlobals::reduceTime_ << std::endl
//     << "\twait     : " << PstreamGlobals::waitTime_ << std::endl
//     << "\tgather   : " << PstreamGlobals::gatherTime_ << std::endl
//     << "\tscatter  : " << PstreamGlobals::scatterTime_ << std::endl
//     << "\tallToAll : " << PstreamGlobals::allToAllTime_ << std::endl
//     << std::endl;


    return true;
}


bool Foam::functionObjects::parProfiling::write()
{
    return true;
}


// ************************************************************************* //
