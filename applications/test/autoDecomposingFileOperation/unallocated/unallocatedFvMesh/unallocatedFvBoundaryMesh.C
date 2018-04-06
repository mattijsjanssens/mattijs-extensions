/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "unallocatedFvBoundaryMesh.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::unallocatedFvBoundaryMesh::findPatchID
(
    const word& patchName
) const
{
    const PtrList<unallocatedGenericFvPatch>& patches = *this;

    forAll(patches, patchi)
    {
        if (patches[patchi].name() == patchName)
        {
            return patchi;
        }
    }

    // Not found, return -1
    return -1;
}


Foam::label Foam::unallocatedFvBoundaryMesh::whichPatch
(
    const label faceIndex
) const
{
    const PtrList<unallocatedGenericFvPatch>& patches = *this;

    forAll(patches, patchi)
    {
        const unallocatedGenericFvPatch& pp = patches[patchi];

        if (pp.start() <= faceIndex && (pp.start()+pp.size() > faceIndex))
        {
            return patchi;
        }
    }

    // Not found, return -1
    return -1;
}


Foam::labelList Foam::unallocatedFvBoundaryMesh::findIndices
(
    const keyType& key,
    const bool usePatchGroups
) const
{
    DynamicList<label> indices;

    if (!key.empty())
    {
        const PtrList<unallocatedGenericFvPatch>& patches = *this;

        wordList names(patches.size());
        forAll(patches, patchi)
        {
            names[patchi] = patches[patchi].name();
        }

        if (key.isPattern())
        {
            indices = findStrings(key, names);

            // if (usePatchGroups && groupPatchIDs().size())
            // {
            //     labelHashSet indexSet(indices);
            //
            //     const wordList allGroupNames = groupPatchIDs().toc();
            //     labelList groupIDs = findStrings(key, allGroupNames);
            //     forAll(groupIDs, i)
            //     {
            //         const word& grpName = allGroupNames[groupIDs[i]];
            //         const labelList& patchIDs = groupPatchIDs()[grpName];
            //         forAll(patchIDs, j)
            //         {
            //             if (indexSet.insert(patchIDs[j]))
            //             {
            //                 indices.append(patchIDs[j]);
            //             }
            //         }
            //     }
            // }
        }
        else
        {
            // Literal string. Special version of above to avoid
            // unnecessary memory allocations

            indices.setCapacity(1);
            forAll(patches, i)
            {
                if (key == patches[i].name())
                {
                    indices.append(i);
                    break;
                }
            }

            //if (usePatchGroups && groupPatchIDs().size())
            //{
            //    const HashTable<labelList, word>::const_iterator iter =
            //        groupPatchIDs().find(key);
            //
            //    if (iter != groupPatchIDs().end())
            //    {
            //        labelHashSet indexSet(indices);
            //
            //        const labelList& patchIDs = iter();
            //        forAll(patchIDs, j)
            //        {
            //            if (indexSet.insert(patchIDs[j]))
            //            {
            //                indices.append(patchIDs[j]);
            //            }
            //        }
            //    }
            //}
        }
    }

    return indices;
}


// ************************************************************************* //
