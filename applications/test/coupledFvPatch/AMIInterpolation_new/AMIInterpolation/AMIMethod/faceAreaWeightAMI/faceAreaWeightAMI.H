/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::faceAreaWeightAMI

Description
    Face area weighted Arbitrary Mesh Interface (AMI) method

SourceFiles
    faceAreaWeightAMI.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaWeightAMI_H
#define faceAreaWeightAMI_H

#include "AMIMethod.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class faceAreaWeightAMI Declaration
\*---------------------------------------------------------------------------*/

class faceAreaWeightAMI
:
    public AMIMethod
{
private:

    // Private data

        //- Flag to restart uncovered source faces
        const bool restartUncoveredSourceFace_;


protected:

    // Protected Member Functions

        //- Disallow default bitwise copy construct
        faceAreaWeightAMI(const faceAreaWeightAMI&);

        //- Disallow default bitwise assignment
        void operator=(const faceAreaWeightAMI&);

        // Marching front

            //- Calculate addressing and weights using temporary storage
            virtual void calcAddressing
            (
                List<DynamicList<label>>& srcAddress,
                List<DynamicList<scalar>>& srcWeights,
                List<DynamicList<point>>& srcCentroid,
                List<DynamicList<label>>& tgtAddress,
                List<DynamicList<scalar>>& tgtWeights,
                List<DynamicList<point>>& tgtCentroid,
                label srcFacei,
                label tgtFacei
            );

            //- Determine overlap contributions for source face srcFacei
            virtual bool processSourceFace
            (
                const label srcFacei,
                const label tgtStartFacei,
                DynamicList<label>& nbrFaces,
                DynamicList<label>& visitedFaces,
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<DynamicList<point>>& srcCentroid,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght,
                List<DynamicList<point>>& tgtCentroid
            );

            //- Attempt to re-evaluate source faces that have not been included
            virtual void restartUncoveredSourceFace
            (
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<DynamicList<point>>& srcCentroid,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght,
                List<DynamicList<point>>& tgtCentroid
            );

            //- Set the source and target seed faces
            virtual void setNextFaces
            (
                label& startSeedI,
                label& srcFacei,
                label& tgtFacei,
                const boolList& mapFlag,
                labelList& seedFaces,
                const DynamicList<label>& visitedFaces,
                bool errorOnNotFound = true
            ) const;


        // Evaluation

            //- The minimum weight below which connections are discarded
            virtual scalar minWeight() const;

            //- Area of intersection between source and target faces
            virtual Tuple2<scalar, point> interArea
            (
                const label srcFacei,
                const label tgtFacei
            ) const;


public:

    //- Runtime type information
    TypeName("faceAreaWeightAMI");


    // Constructors

        //- Construct from components
        faceAreaWeightAMI
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true,
            const bool restartUncoveredSourceFace = true
        );


    //- Destructor
    virtual ~faceAreaWeightAMI();


    // Member Functions

        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                List<List<point>>& srcCentroids,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                List<List<point>>& tgtCentroids,
                label srcFacei = -1,
                label tgtFacei = -1
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
