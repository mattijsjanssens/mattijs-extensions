/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::cyclicAMILduInterface

Description
    An abstract base class for cyclic AMI coupled interfaces

SourceFiles
    cyclicAMILduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMILduInterface_H
#define cyclicAMILduInterface_H

#include "primitiveFieldsFwd.H"
#include "AMIPatchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicAMILduInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMILduInterface
{

public:

    //- Runtime type information
    TypeName("cyclicAMILduInterface");


    // Constructors

        //- Construct null
        cyclicAMILduInterface() = default;


    //- Destructor
    virtual ~cyclicAMILduInterface() = default;


    // Member Functions

        // Access

            //- Neighbour patch ID
            virtual const labelList& neighbPatchIDs() const = 0;

            //- Am I owner of this 'connection'
            virtual bool owner(const label index) const = 0;

            //- Return neighbour
            virtual const cyclicAMILduInterface& neighbPatch
            (
                const label index
            ) const = 0;

            //- Return AMI (or null)
            virtual tmpNrc<AMIPatchToPatchInterpolation> AMI
            (
                const label index
            ) const = 0;

            //- Total neighbour size (of neighbours with valid AMI only)
            virtual label neighbSize() const = 0;

            //- Indices of neighbours with valid AMI
            virtual const labelList& AMIIndices() const = 0;

            //- Index on neighbour interface of my interface
            virtual label neighbIndex(const label index) const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return face reverse transformation tensor
            virtual const tensorField& reverseT() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
