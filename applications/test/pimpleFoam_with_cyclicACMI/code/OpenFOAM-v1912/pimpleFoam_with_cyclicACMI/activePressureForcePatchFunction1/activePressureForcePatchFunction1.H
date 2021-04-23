/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 M. Janssens
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
    Foam::PatchFunction1Types::activePressureForce

Description
    See activePressureForceBaffleVelocityFvPatchVectorField but in
    PatchFunction1 for use in a cyclicACMI form.

SourceFiles
    activePressureForcePatchFunction1.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFunction1Types_activePressureForce_H
#define PatchFunction1Types_activePressureForce_H

#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                       Class activePressureForce Declaration
\*---------------------------------------------------------------------------*/

class activePressureForce
:
    public PatchFunction1<scalar>
{
    // Private data

        //- Name of the pressure field used to calculate the force
        //- on the active baffle
        word pName_;

//        //- Name of the cyclic patch used when the active baffle is open
//        word cyclicPatchName_;
//
//        //- Index of the cyclic patch used when the active baffle is open
//        label cyclicPatchLabel_;

        //- Initial wall patch areas
        vectorField initWallSf_;

        //- Initial cyclic patch areas
        vectorField initCyclicSf_;

        //- Initial neighbour-side cyclic patch areas
        vectorField nbrCyclicSf_;

        //- Current fraction of the active baffle which is open
        scalar openFraction_;

        //- Time taken for the active baffle to open
        scalar openingTime_;

        //- Maximum fractional change to the active baffle openness
        //- per time-step
        scalar maxOpenFractionDelta_;

        label curTimeIndex_;

        //- Minimum value for the active baffle to start opening
        scalar minThresholdValue_;

        //- Force based active baffle
        bool fBased_;

        //- Baffle is activated
        bool baffleActivated_;

        //- Baffle is opening or closing (1 opening, 0 closing)
        bool opening_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const activePressureForce&) = delete;


public:

    // Runtime type information
    TypeName("activePressureForce");


    // Constructors

        //- Construct from entry name and dictionary
        activePressureForce
        (
            const polyPatch& pp,
            const word& type,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy constructor
        explicit activePressureForce(const activePressureForce& ut);

        //- Copy constructor setting patch
        explicit activePressureForce
        (
            const activePressureForce& ut,
            const polyPatch& pp
        );

        //- Construct and return a clone
        virtual tmp<PatchFunction1> clone() const
        {
            return tmp<PatchFunction1<scalar>>
            (
                new activePressureForce(*this)
            );
        }

        //- Construct and return a clone setting patch
        virtual tmp<PatchFunction1<scalar>> clone(const polyPatch& pp) const
        {
            return tmp<PatchFunction1<scalar>>
            (
                new activePressureForce(*this, pp)
            );
        }


    //- Destructor
    virtual ~activePressureForce() = default;


    // Member Functions

        // Evaluation

            //- Return activePressureForce value
            virtual tmp<Field<scalar>> value(const scalar x) const;

            //- Is value constant (i.e. independent of x)
            virtual bool constant() const
            {
                return false;
            }

            //- Is value uniform (i.e. independent of coordinate)
            virtual bool uniform() const
            {
                return true;
            }

            //- Integrate between two values
            virtual tmp<Field<scalar>> integrate
            (
                const scalar x1,
                const scalar x2
            ) const
            {
                notImplemented;
                return tmpField<scalar>>::null;
            }


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const FieldMapper& mapper);

            //- Reverse map the given PatchFunction1 onto this PatchFunction1
            virtual void rmap
            (
                const PatchFunction1<scalar>& pf1,
                const labelList& addr
            );


        // I-O

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "activePressureForceI.H"
//
//#ifdef NoRepository
//    #include "activePressureForce.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
