/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fv::gaussDivScheme2

Group
    grpFvDivSchemes

Description
    Basic second-order div using face-gradients and Gauss' theorem.

SourceFiles
    gaussDivScheme2.C

\*---------------------------------------------------------------------------*/

#ifndef gaussDivScheme2_H
#define gaussDivScheme2_H

#include "divScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class gaussDivScheme2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussDivScheme2
:
    public fv::divScheme<Type>
{
    // Private Member Functions

        //- No copy construct
        gaussDivScheme2(const gaussDivScheme2&) = delete;

        //- No copy assignment
        void operator=(const gaussDivScheme2&) = delete;


public:

    //- Runtime type information
    TypeName("Gauss2");


    // Constructors

        //- Construct null
        gaussDivScheme2(const fvMesh& mesh)
        :
            divScheme<Type>(mesh)
        {}

        //- Construct from Istream
        gaussDivScheme2(const fvMesh& mesh, Istream& is)
        :
            divScheme<Type>(mesh, is)
        {}


    // Member Functions

        tmp
        <
            GeometricField
            <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
        > fvcDiv
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

/*
        //- Inplace
        void fvcDiv
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& lambdas,
            GeometricField
            <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>&
        );
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gaussDivScheme2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
