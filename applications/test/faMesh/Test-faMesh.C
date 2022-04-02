/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Group

Description

Author

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Passive scalar transport finiteArea equation solver."
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFaMesh.H"

    const auto& aPatch = aMesh.patch();
    Pout<< "nFaces:" << aPatch.size()
        << " nInternalEdges:" << aPatch.nInternalEdges() << endl;


    // faceAreaNormals  : normals of boundary faces
    // edgeAreaNormals  : average of two neighbouring boundary faces
    // Le               : outer product of edgeVector and edgeAreaNormal

    vectorField faceAreaNormals(mesh.faceAreas(), aMesh.faceLabels());
    faceAreaNormals /= mag(faceAreaNormals);
    vectorField edgeAreaNormals;
    vectorField Le;
    {
        const auto& edgeFaces = aPatch.edgeFaces();
        edgeAreaNormals.setSize(edgeFaces.size(), Zero);
        forAll(edgeFaces, edgei)
        {
            const auto& eFaces = edgeFaces[edgei];
            for (const label facei : eFaces)
            {
                edgeAreaNormals[edgei] += faceAreaNormals[facei];
            }
            edgeAreaNormals[edgei] /= mag(edgeAreaNormals[edgei]);
        }

        Le.setSize(edgeAreaNormals.size());
        forAll(Le, edgei)
        {
            const edge& e = aPatch.edges()[edgei];
            Le[edgei] = e.vec(aPatch.localPoints()) ^ edgeAreaNormals[edgei];
        }
    }





    Pout<< "Mesh:" << nl
        << incrIndent
        << indent << "nFaces:" << aMesh.nFaces() << nl
        << indent << "Le:" << aMesh.Le() << nl
        << indent << "MY Le:" << Le << nl
        << indent << "magLe:" << aMesh.magLe() << nl
        << indent << "areaCentres:" << aMesh.areaCentres() << nl
        << indent << "faceAreas:" << aPatch.faceCentres() << nl
        << indent << "edgeCentres:" << aMesh.edgeCentres() << nl
        << indent << "S:" << aMesh.S() << nl
        << indent << "faceAreaNormals:" << aMesh.faceAreaNormals() << nl
        << indent << "MY faceAreaNormals:" << faceAreaNormals << nl
        << indent << "edgeAreaNormals:" << aMesh.edgeAreaNormals() << nl
        << indent << "MY edgeAreaNormals:" << edgeAreaNormals << nl
        << indent << "pointAreaNormals:" << aMesh.pointAreaNormals() << nl
        << indent << "faceCurvatures:" << aMesh.faceCurvatures() << nl
        << decrIndent
        << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
