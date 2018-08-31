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
    Foam::unallocatedProcessorFvPatchField

Description

SourceFiles
    unallocatedProcessorFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef unallocatedProcessorFvPatchField_H
#define unallocatedProcessorFvPatchField_H

#include "unallocatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class unallocatedProcessorFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class unallocatedProcessorFvPatchField
:
    public unallocatedFvPatchField<Type>
{
public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from patch and internal field
        unallocatedProcessorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unallocatedProcessorFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        unallocatedProcessorFvPatchField
        (
            const unallocatedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, unallocatedVolMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        unallocatedProcessorFvPatchField
        (
            const unallocatedProcessorFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<unallocatedFvPatchField<Type>> clone() const
        {
            return tmp<unallocatedFvPatchField<Type>>
            (
                new unallocatedProcessorFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unallocatedProcessorFvPatchField
        (
            const unallocatedProcessorFvPatchField<Type>&,
            const DimensionedField<Type, unallocatedVolMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<unallocatedFvPatchField<Type>> clone
        (
            const DimensionedField<Type, unallocatedVolMesh>& iF
        ) const
        {
            return tmp<unallocatedFvPatchField<Type>>
            (
                new unallocatedProcessorFvPatchField<Type>(*this, iF)
            );
        }


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unallocatedProcessorFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
