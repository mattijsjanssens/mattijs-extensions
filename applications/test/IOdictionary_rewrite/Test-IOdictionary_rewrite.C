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

Application

Description

\*---------------------------------------------------------------------------*/

#include "Field.H"
#include "fvMesh.H"
#include "argList.H"
#include "fieldDictionary.H"
#include "volFields.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void replaceToken(dictionary& dict, const token& oldTok, const token& newTok)
{
    forAllIter(IDLList<entry>, dict, iter)
    {
        if (iter().isDict())
        {
            Pout<< incrIndent;
            replaceToken(iter().dict(), oldTok, newTok);
            Pout<< decrIndent;
        }
        else
        {
            ITstream& str = iter().stream();
            forAll(str, i)
            {
                Pout<< indent << "i:" << i << endl;
                Pout<< indent << str[i].info() << endl;
                Pout<< endl;

                if (str[i] == oldTok)
                {
                    str[i] = newTok;
                }
            }
        }
    }
}



int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Try to load field
    IOobject fieldObject
    (
        "p",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
    );

    DebugVar(fieldObject.headerClassName());

    IFstream str(typeFilePath<labelIOList>(fieldObject));

    // Read dictionary
    fieldDictionary fieldDict
    (
        fieldObject,
        "volScalarField"
    );

    DebugVar(fieldDict);

    replaceToken
    (
        fieldDict,
        token(word("nonuniform")),
        token(word("unchecked"))
    );

    DebugVar(fieldDict);

    Field<scalar> fld("internalField", fieldDict, 1);
    DebugVar(fld);

    Info<< "end" << endl;
}


// ************************************************************************* //
