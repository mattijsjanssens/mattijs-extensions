/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "fvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::pTraits<Type>::labelType Foam::fvMesh::validComponents() const
{
    return pow
    (
        this->solutionD(),
        pTraits
        <
            typename powProduct<Vector<label>,
            pTraits<Type>::rank>::type
        >::zero
    );
}


template<class GeoField>
void Foam::fvMesh::addPatchFields
(
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const typename GeoField::value_type& defaultPatchValue
)
{
    HashTable<GeoField*> flds(thisDb().lookupClass<GeoField>());

    const wordList fldNames(flds.sortedToc());

    forAll(fldNames, i)
    {
        GeoField& fld = *flds[fldNames[i]];

        typename GeoField::Boundary& bfld =
            fld.boundaryFieldRef();

        label sz = bfld.size();
        bfld.setSize(sz+1);

        if (patchFieldDict.found(fld.name()))
        {
            bfld.set
            (
                sz,
                GeoField::Patch::New
                (
                    fld.mesh().boundary()[sz],
                    fld(),
                    patchFieldDict.subDict(fld.name())
                )
            );
        }
        else
        {
            bfld.set
            (
                sz,
                GeoField::Patch::New
                (
                    defaultPatchFieldType,
                    fld.mesh().boundary()[sz],
                    fld()
                )
            );
            bfld[sz] == defaultPatchValue;
        }
    }
}


//template<class GeoField>
//void Foam::fvMesh::setPatchFields
//(
//    const label patchi,
//    const dictionary& patchFieldDict
//)
//{
//    HashTable<GeoField*> flds(thisDb().lookupClass<GeoField>());
//
//    const wordList fldNames(flds.sortedToc());
//
//    forAll(fldNames, i)
//    {
//        GeoField& fld = *flds[fldNames[i]];
//
//        typename GeoField::Boundary& bfld =
//            fld.boundaryFieldRef();
//
//        if (patchFieldDict.found(fld.name()))
//        {
//            bfld.set
//            (
//                patchi,
//                GeoField::Patch::New
//                (
//                    fld.mesh().boundary()[patchi],
//                    fld(),
//                    patchFieldDict.subDict(fld.name())
//                )
//            );
//        }
//    }
//}
//
//
//template<class GeoField>
//void Foam::fvMeshTools::setPatchFields
//(
//    objectRegistry& obr,
//    const label patchi,
//    const typename GeoField::value_type& value
//)
//{
//    HashTable<GeoField*> flds(obr.lookupClass<GeoField>());
//
//    const wordList fldNames(flds.sortedToc());
//
//    forAll(fldNames, i)
//    {
//        GeoField& fld = *flds[fldNames[i]];
//        fld.boundaryFieldRef()[patchi] == value;
//    }
//}


// Remove last patch field
template<class GeoField>
void Foam::fvMesh::trimPatchFields
(
    const label nPatches
)
{
    HashTable<GeoField*> flds(thisDb().lookupClass<GeoField>());

    const wordList fldNames(flds.sortedToc());

    forAll(fldNames, i)
    {
        GeoField& fld = *flds[fldNames[i]];
        fld.boundaryFieldRef().setSize(nPatches);
    }
}


// Reorder patch field
template<class GeoField>
void Foam::fvMeshTools::reorderPatchFields
(
    const labelList& oldToNew
)
{
    HashTable<GeoField*> flds(thisDb().lookupClass<GeoField>());

    const wordList fldNames(flds.sortedToc());

    forAll(fldNames, i)
    {
        GeoField& fld = *flds[fldNames[i]];
        fld.boundaryFieldRef().reorder(oldToNew);
    }
}


// ************************************************************************* //
