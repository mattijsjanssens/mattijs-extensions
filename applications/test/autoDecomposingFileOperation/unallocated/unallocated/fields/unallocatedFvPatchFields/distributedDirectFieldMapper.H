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

Class
    Foam::distributedDirectFieldMapper

Description
    FieldMapper with direct mapping from remote quantities.

    Use distributedUnallocatedDirectFieldMapper if map does not need to be
    stored.

\*---------------------------------------------------------------------------*/

#ifndef distributedDirectFieldMapper_H
#define distributedDirectFieldMapper_H

#include "Field.H"
#include "FieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class distributedDirectFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class distributedDirectFieldMapper
:
    public FieldMapper
{
    const labelUList& directAddressing_;

    mapDistributeBase distMap_;

    bool hasUnmapped_;

public:

    // Constructors

        //- Construct given addressing. If constructSize is -1 it will
        //  use the distMap as is, otherwise will construct the reverse map.
        distributedDirectFieldMapper
        (
            const labelUList& directAddressing,
            const mapDistributeBase& distMap,
            const label constructSize
        );
//
//        :
//            directAddressing_(directAddressing),
//            hasUnmapped_(false)
//        {
//            if
//            (
//                notNull(directAddressing_)
//             && directAddressing_.size()
//             && min(directAddressing_) < 0
//            )
//            {
//                hasUnmapped_ = true;
//            }
//        }

    //- Destructor
    virtual ~distributedDirectFieldMapper()
    {}


    // Member Functions

        virtual label size() const
        {
            return
            (
                notNull(directAddressing_)
              ? directAddressing_.size()
              : distMap_.constructSize()
            );
        }

        virtual bool direct() const
        {
            return true;
        }

        virtual bool distributed() const
        {
            return true;
        }

        virtual const mapDistributeBase& distributeMap() const
        {
            return distMap_;
        }

        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        virtual const labelUList& directAddressing() const
        {
            return directAddressing_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
