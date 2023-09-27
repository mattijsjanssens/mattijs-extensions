/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 M. Janssens
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
    Test-primitivePatch_distribute

Description
    Showcase distributing primitivePatch

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "OBJstream.H"
#include "DynamicField.H"
#include "cyclicAMIPolyPatch.H"
#include "columnFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Equivalence of mapDistributeBase::distribute
void distribute
(
    const bool forward,
    const mapDistributeBase& distMap,
    const faceList& faces,
    const pointField& points,

    DynamicList<face>& newFaces,
    DynamicField<point>& newPoints
)
{
    // Use resulting map to send over faces/points. Easiest with PstreamBuffers
    // since unpacked data
    PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

    {
        const labelListList& map =
            (forward ? distMap.subMap() : distMap.constructMap());
        forAll(map, proci)
        {
            if (map[proci].size())
            {
                faceList subFaces(faces, map[proci]);
                const primitivePatch subPatch
                (
                    SubList<face>(subFaces, subFaces.size()),
                    points
                );
                UOPstream os(proci, pBufs);

                // TBD: send as compactListList?
                os << subPatch.localFaces() << subPatch.localPoints();
            }
        }
    }
    pBufs.finishedSends();

    newFaces.setCapacity(distMap.constructSize());  // guesstimate
    newPoints.setCapacity(distMap.constructSize()); // guesstimate

    {
        const labelListList& map =
            (forward ? distMap.constructMap() : distMap.subMap());
        forAll(map, proci)
        {
            if (map[proci].size())
            {
                UIPstream is(proci, pBufs);

                // TBD. map not used. Assumed sending in processor order
                faceList subFaces;
                pointField subPoints;
                is >> subFaces >> subPoints;

                if (subFaces.size() != map[proci].size())
                {
                    FatalErrorInFunction<< "problem" << exit(FatalError);
                }

                // Append (renumbered) faces
                const label offset = newPoints.size();
                for (auto& f : subFaces)
                {
                    f += offset;
                    newFaces.append(f);
                }
                newPoints.append(subPoints);
            }
        }
    }
}




// Main program:

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const label patchi = 0;

    const auto& pbm = mesh.boundaryMesh();

    const primitivePatch& pp = pbm[patchi];

    const globalIndex globalFaces(pp.size());


    // Decide which face goes where
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList destProc(pp.size());
    forAll(pp, i)
    {
        destProc[i] = UPstream::nProcs()-1-UPstream::myProcNo();
    }


    // Construct map
    // ~~~~~~~~~~~~~

    const mapDistributeBase map(invertOneToMany(UPstream::nProcs(), destProc));

    //Pout<< "subMap:" << map.subMap() << endl;
    //Pout<< "constructMap:" << map.constructMap() << endl;


    // Send over faces+points
    // ~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<face> newFaces;
    DynamicField<point> newPoints;
    distribute
    (
        true,       // send
        map,
        pp,
        pp.points(),

        newFaces,
        newPoints
    );

    //Pout<< "newFaces:" << flatOutput(newFaces) << endl;
    //Pout<< "newPoints:" << flatOutput(newPoints) << endl;


    // Do something with newFaces,newPoints
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointField newCentres;
    {
        const primitivePatch newPatch
        (
            SubList<face>(newFaces, newFaces.size()),
            newPoints
        );
        newCentres = newPatch.faceCentres();

        OBJstream os(runTime.path()/"newPatch.obj");
        os.write(newFaces, newPoints, false);
        Pout<< "Written " << newFaces.size() << " faces, " << newPoints.size()
            << " points to " << os.name() << endl;
    }


    // Send back
    // ~~~~~~~~~

    map.reverseDistribute(pp.size(), newCentres);

    forAll(newCentres, i)
    {
        Pout<< "face:" << i << " old:" << pp.faceCentres()[i]
            << " new:" << newCentres[i]
            << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
