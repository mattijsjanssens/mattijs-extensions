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

\*---------------------------------------------------------------------------*/

#ifndef unallocatedProcessorFvsPatchFields_H
#define unallocatedProcessorFvsPatchFields_H

#include "unallocatedProcessorFvsPatchField.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef unallocatedProcessorFvsPatchField<scalar>
    processorUnallocatedFvsPatchScalarField;
typedef unallocatedProcessorFvsPatchField<vector>
    processorUnallocatedFvsPatchVectorField;
typedef unallocatedProcessorFvsPatchField<sphericalTensor>
    processorUnallocatedFvsPatchSphericalTensorField;
typedef unallocatedProcessorFvsPatchField<symmTensor>
    processorUnallocatedFvsPatchSymmTensorField;
typedef unallocatedProcessorFvsPatchField<tensor>
    processorUnallocatedFvsPatchTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
