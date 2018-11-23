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

\*---------------------------------------------------------------------------*/

#include "subsetTriSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "triSurfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(subsetTriSurfaceMesh, 0);
    addToRunTimeSelectionTable(searchableSurface, subsetTriSurfaceMesh, dict);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Foam::fileName Foam::subsetTriSurfaceMesh::checkFile
// (
//     const regIOobject& io,
//     const bool isGlobal
// )
// {
//     const fileName fName
//     (
//         isGlobal
//       ? io.globalFilePath(typeName)
//       : io.localFilePath(typeName)
//     );
//     if (fName.empty())
//     {
//         FatalErrorInFunction
//             << "Cannot find subsetTriSurfaceMesh starting from "
//             << io.objectPath() << exit(FatalError);
//     }
// 
//     return fName;
// }
// 
// 
// Foam::fileName Foam::subsetTriSurfaceMesh::relativeFilePath
// (
//     const regIOobject& io,
//     const fileName& f,
//     const bool isGlobal
// )
// {
//     fileName fName(f);
//     fName.expand();
//     if (!fName.isAbsolute())
//     {
//         // Is the specified file:
//         // - local to the cwd?
//         // - local to the case dir?
//         // - or just another name?
//         fName = fileHandler().filePath
//         (
//             isGlobal,
//             IOobject(io, fName),
//             word::null
//         );
//     }
//     return fName;
// }
// 
// 
// Foam::fileName Foam::subsetTriSurfaceMesh::checkFile
// (
//     const regIOobject& io,
//     const dictionary& dict,
//     const bool isGlobal
// )
// {
//     fileName dictFName, fName;
// 
//     if (dict.readIfPresent("file", dictFName, false, false))
//     {
//         fName = relativeFilePath(io, dictFName, isGlobal);
// 
//         if (!exists(fName))
//         {
//             FatalErrorInFunction
//                 << "Cannot find subsetTriSurfaceMesh at " << io.path(dictFName)
//                 << exit(FatalError);
//         }
//     }
//     else
//     {
//         fName =
//         (
//             isGlobal
//           ? io.globalFilePath(typeName)
//           : io.localFilePath(typeName)
//         );
// 
//         if (!exists(fName))
//         {
//             FatalErrorInFunction
//                 << "Cannot find subsetTriSurfaceMesh starting from "
//                 << io.objectPath() << exit(FatalError);
//         }
//     }
// 
//     return fName;
// }

Foam::triSurface Foam::subsetTriSurfaceMesh::subsetTriSurface
(
    const IOobject& io,
    const dictionary& dict
)
{
DebugVar(dict);
    //fileName fName(io.globalFilePath(triSurfaceMesh::typeName));
    fileName fName(triSurfaceMesh::checkFile(io, dict, true));

DebugVar(fName);

    // Read surface
    triSurface s(fName);

    const word fieldName(dict.lookup("field"));

    // Read field
    triSurfaceScalarField fld
    (
        IOobject
        (
            fieldName,
            io.time().constant(),
            "triSurface",
            io.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        s
    );

    // Subset mesh
    const scalar min(readScalar(dict.lookup("min")));
    const scalar max(readScalar(dict.lookup("max")));

    boolList include(s.size());
    forAll(fld, i)
    {
        include[i] = (fld[i] >= min && fld[i] <= max);
    }

//DebugVar(include);

    labelList pointMap, faceMap;
    triSurface sub(s.subsetMesh(include, pointMap, faceMap));
    sub.writeStats(Info);

    return sub;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subsetTriSurfaceMesh::subsetTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, subsetTriSurface(io, dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::subsetTriSurfaceMesh::~subsetTriSurfaceMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
