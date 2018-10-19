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

Class
    Foam::basicParticle

Description
    Storage and IO for particle without mesh

\*---------------------------------------------------------------------------*/

#ifndef basicParticle_H
#define basicParticle_H

#include "vector.H"
#include "barycentric.H"
#include "IDLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class basicParticle;

// Forward declaration of friend functions and operators

Ostream& operator<<
(
    Ostream&,
    const basicParticle&
);

/*---------------------------------------------------------------------------*\
                        Class basicParticle Declaration
\*---------------------------------------------------------------------------*/

class basicParticle
:
    public IDLList<basicParticle>::link
{
    // Private data

        //- Size in bytes of the position data
        static const std::size_t sizeofPosition_;

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- Coordinates of particle
        barycentric coordinates_;

        //- Index of the cell it is in
        label celli_;

        //- Index of the face that owns the decomposed tet that the
        //  particle is in
        label tetFacei_;

        //- Index of the point on the face that defines the decomposed
        //  tet that the particle is in.  Relative to the face base
        //  point.
        label tetPti_;

        //- Face index if the particle is on a face otherwise -1
        label facei_;

        //- Fraction of time-step completed
        scalar stepFraction_;

        //- Originating processor id
        label origProc_;

        //- Local particle id on originating processor
        label origId_;


public:

    // Static data members

         //- Runtime type information
         ClassName("particle");


    // Constructors

        //- Construct from Istream
        basicParticle(Istream&, bool readFields = true);

        //- Construct copy with new addressing
        basicParticle
        (
            const basicParticle&,
            const label celli,
            const label tetFacei,
            const label tetPti
        );


    //- Destructor
    virtual ~basicParticle()
    {}


    // Member Functions

        // Access

            //- Return current particle coordinates
            inline const barycentric& coordinates() const;

            //- Return current cell particle is in
            inline label cell() const;

            //- Return current tet face particle is in
            inline label tetFace() const;

            //- Return current tet face particle is in
            inline label tetPt() const;

            //- Return current face particle is on otherwise -1
            inline label face() const;

            //- Return the fraction of time-step completed
            inline scalar stepFraction() const;

            //- Return the originating processor ID
            inline label origProc() const;

            //- Return the particle ID on the originating processor
            inline label origId() const;


    // I-O

        static autoPtr<basicParticle> New(Istream& is)
        {
            return autoPtr<basicParticle>(new basicParticle(is, false));
        }

        //- Write the particle position and cell
        void writePosition(Ostream&) const;


    // Friend Operators

        friend Ostream& operator<<(Ostream&, const basicParticle&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
