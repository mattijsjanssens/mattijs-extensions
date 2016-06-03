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
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "SHA1Digest.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
//#include "stringOps.H"
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
    dynCode.addCompileFile("functionObjectTemplate.C");

    // Copy filtered H template
    dynCode.addCopyFile("functionObjectTemplate.H");

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
      + context.filterVars()["codeOptions"]
      + "\n\nLIB_LIBS = \\\n"
      + "    -lOpenFOAM \\\n"
      + "    -lfiniteVolume \\\n"
      + "    -lmeshTools \\\n"
      + context.filterVars()["codeLibs"]
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


bool Foam::codedFunctionObject::execute(const bool postProcess)
{
    updateLibrary(name_, context_);
    return redirectFunctionObject().execute(postProcess);
}


bool Foam::codedFunctionObject::write(const bool postProcess)
{
    updateLibrary(name_, context_);
    return redirectFunctionObject().write(postProcess);
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


    // Extract all filter variables from dictionary and calculate sha1.
    context_.clear();

    //- Note: could assume all strings in dictionary are filter variables
    //        but this would also include e.g. the write settings
    //context_.read(dict);

    // From looking through functionObjectTemplate.[CH] :

    context_.addFilterVariable(false, dict, "codeOptions");
    context_.addFilterVariable(false, dict, "codeLibs");

    context_.addFilterVariable(true, dict, "codeInclude");
    context_.addFilterVariable(true, dict, "codeData");
    context_.addFilterVariable(true, dict, "localCode");
    context_.addFilterVariable(true, dict, "codeRead");
    context_.addFilterVariable(true, dict, "codeExecute");
    context_.addFilterVariable(true, dict, "codeWrite");
    context_.addFilterVariable(true, dict, "codeEnd");

    updateLibrary(name_, context_);
    return redirectFunctionObject().read(dict);
}


// ************************************************************************* //
