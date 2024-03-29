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

    //- Reduction class. If x and y are not equal assign value.
    class statsEqOp
    {
        public:
        void operator()
        (
            FixedList<statData, 2>& xStats,
            const FixedList<statData, 2>& yStats
        ) const
        {
            forAll(xStats, i)
            {
                statData& x = xStats[i];
                const statData& y = yStats[i];

                // 0 : min
                // 1 : max
                // 2 : sum
                if (y[0].second() < x[0].second())
                {
                    x[0].second() = y[0].second();
                    x[0].first()  = y[0].first();
                }
                if (y[1].second() > x[1].second())
                {
                    x[1].second() = y[1].second();
                    x[1].first()  = y[1].first();
                }
                x[2].second() += y[2].second();
                x[2].first()++;
            }
        }
    };
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
    typedef FixedList<Tuple2<label, scalar>, 3> statData;
    FixedList<statData, 2> times;

    {
        const scalar masterTime =
            PstreamGlobals::reduceTime_
          + PstreamGlobals::gatherTime_
          + PstreamGlobals::scatterTime_;

        statData& reduceStats = times[0];

        Tuple2<label, scalar>& minTime = reduceStats[0];
        minTime.first() = Pstream::myProcNo();
        minTime.second() = masterTime;

        Tuple2<label, scalar>& maxTime = reduceStats[1];
        maxTime.first() = Pstream::myProcNo();
        maxTime.second() = masterTime;

        Tuple2<label, scalar>& sumTime = reduceStats[2];
        sumTime.first() = 1;
        sumTime.second() = masterTime;
    }
    {
        const scalar allTime =
            PstreamGlobals::waitTime_
          + PstreamGlobals::allToAllTime_;

        statData& allToAllStats = times[1];

        Tuple2<label, scalar>& minTime = allToAllStats[0];
        minTime.first() = Pstream::myProcNo();
        minTime.second() = allTime;

        Tuple2<label, scalar>& maxTime = allToAllStats[1];
        maxTime.first() = Pstream::myProcNo();
        maxTime.second() = allTime;

        Tuple2<label, scalar>& sumTime = allToAllStats[2];
        sumTime.first() = 1;
        sumTime.second() = allTime;
    }


    {
        const scalar oldReduceTime(PstreamGlobals::reduceTime_);
        const scalar oldWaitTime(PstreamGlobals::waitTime_);
        const scalar oldGatherTime(PstreamGlobals::gatherTime_);
        const scalar oldScatterTime(PstreamGlobals::scatterTime_);
        const scalar oldAllToAlltime(PstreamGlobals::allToAllTime_);

        Pstream::combineGather(times, statsEqOp());

        PstreamGlobals::reduceTime_ = oldReduceTime;
        PstreamGlobals::waitTime_ = oldWaitTime;
        PstreamGlobals::gatherTime_ = oldGatherTime;
        PstreamGlobals::scatterTime_ = oldScatterTime;
        PstreamGlobals::allToAllTime_ = oldAllToAlltime;
    }

    if (Pstream::master())
    {
        const statData& reduceStats = times[0];
        const statData& allToAllStats = times[1];

        scalar reduceAvg = reduceStats[2].second()/Pstream::nProcs();
        scalar reduceRelMin = (reduceStats[0].second()-reduceAvg)/reduceAvg;
        scalar reduceRelMax = (reduceStats[1].second()-reduceAvg)/reduceAvg;

        scalar allToAllAvg = allToAllStats[2].second()/Pstream::nProcs();
        scalar allToAllRelMin =
            (allToAllStats[0].second()-allToAllAvg)/allToAllAvg;
        scalar allToAllRelMax =
            (allToAllStats[1].second()-allToAllAvg)/allToAllAvg;

        Info<< type() << ':' << nl
            << "\treduce    : avg = " << reduceAvg << 's' << nl
            << "\t            min = " << reduceRelMin*100
            << "% (on processor " << reduceStats[0].first() << ")" << nl
            << "\t            max = +" << reduceRelMax*100
            << "% (on processor " << reduceStats[1].first() << ")" << nl
            << "\tall-all   : avg = " << allToAllAvg << 's' << nl
            << "\t            min = " << allToAllRelMin*100
            << "% (on processor " << allToAllStats[0].first() << ")" << nl
            << "\t            max = +" << allToAllRelMax*100
            << "% (on processor " << allToAllStats[1].first() << ")" << endl;
    }

    return true;
}


bool Foam::functionObjects::parProfiling::write()
{
    return true;
}


// ************************************************************************* //
