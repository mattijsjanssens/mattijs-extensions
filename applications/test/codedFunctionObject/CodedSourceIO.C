/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "CodedSource.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::CodedSource<Type>::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("fields") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Backward compatibility
        if (dict.found("redirectType"))
        {
            dict.lookup("redirectType") >> name_;
        }
        else
        {
            dict.lookup("name") >> name_;
        }

        // Compilation options
        context_.addFilterVariable(false, dict, "codeOptions");
        context_.addFilterVariable(false, dict, "codeLibs");

        // From looking through the functionObjectTemplate*[CH] :
        context_.addFilterVariables
        (
            dynamicCode::resolveTemplate(codeTemplateC),
            dict
        );
        context_.addFilterVariables
        (
            dynamicCode::resolveTemplate(codeTemplateH),
            dict
        );

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
