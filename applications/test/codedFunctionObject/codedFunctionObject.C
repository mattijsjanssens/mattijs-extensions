/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "codedFunctionObject.H"
#include "Time.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(codedFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        codedFunctionObject,
        dictionary
    );
}

const Foam::word Foam::codedFunctionObject::codeTemplateC
    = "functionObjectTemplate.C";

const Foam::word Foam::codedFunctionObject::codeTemplateH
    = "functionObjectTemplate.H";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::codedFunctionObject::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name_);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    // Debugging: make BC verbose
    // dynCode.setFilterVariable("verbose", "true");
    // Info<<"compile " << name_ << " sha1: "
    //     << context.sha1() << endl;

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
      + dynCode.filterVar("codeOptions")
      + "\n\nLIB_LIBS = \\\n"
      + "    -lOpenFOAM \\\n"
      + "    -lfiniteVolume \\\n"
      + "    -lmeshTools \\\n"
      + dynCode.filterVar("codeLibs")
    );
}


Foam::dlLibraryTable& Foam::codedFunctionObject::libs() const
{
    return const_cast<Time&>(time_).libs();
}


Foam::string Foam::codedFunctionObject::description() const
{
    return "functionObject " + name();
}


void Foam::codedFunctionObject::clearRedirect() const
{
    redirectFunctionObjectPtr_.clear();
}


const Foam::dictionary& Foam::codedFunctionObject::codeDict() const
{
    return dict_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedFunctionObject::codedFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject(name),
    codedBase(),
    time_(time),
    dict_(dict)
{
    read(dict_);

    updateLibrary(name_, context_);
    redirectFunctionObject();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::codedFunctionObject::~codedFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::functionObject& Foam::codedFunctionObject::redirectFunctionObject() const
{
    if (!redirectFunctionObjectPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", name_);

        redirectFunctionObjectPtr_ = functionObject::New
        (
            name_,
            time_,
            constructDict
        );
    }
    return redirectFunctionObjectPtr_();
}


bool Foam::codedFunctionObject::execute()
{
    updateLibrary(name_, context_);
    return redirectFunctionObject().execute();
}


bool Foam::codedFunctionObject::write()
{
    updateLibrary(name_, context_);
    return redirectFunctionObject().write();
}


bool Foam::codedFunctionObject::end()
{
    updateLibrary(name_, context_);
    return redirectFunctionObject().end();
}


bool Foam::codedFunctionObject::read(const dictionary& dict)
{
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

    updateLibrary(name_, context_);
    return redirectFunctionObject().read(dict);
}


// ************************************************************************* //
