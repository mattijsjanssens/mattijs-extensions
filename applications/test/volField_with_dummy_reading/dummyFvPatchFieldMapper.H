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

Class
    Foam::dummyFvPatchFieldMapper

Description
    fvPatchFieldMapper that maps onto zero-sized patch

\*---------------------------------------------------------------------------*/

#ifndef dummyFvPatchFieldMapper_H
#define dummyFvPatchFieldMapper_H

#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class dummyFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class dummyFvPatchFieldMapper
:
    public fvPatchFieldMapper
{
    labelList directAddressing_;

public:

    // Constructors

        //- Construct from components
        dummyFvPatchFieldMapper()
        :
            directAddressing_(0)
        {}


    //- Destructor
    virtual ~dummyFvPatchFieldMapper()
    {}


    // Member Functions

        virtual label size() const
        {
            return directAddressing_.size();
        }

        virtual bool direct() const
        {
            return true;
        }

        virtual bool hasUnmapped() const
        {
            return false;
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
