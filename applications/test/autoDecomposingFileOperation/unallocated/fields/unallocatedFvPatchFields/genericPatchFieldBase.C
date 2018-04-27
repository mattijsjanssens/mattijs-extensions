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

#include "genericPatchFieldBase.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::genericPatchFieldBase::read
(
    const dictionary& dict,
    HashPtrTable<scalarField>& scalarFields,
    HashPtrTable<vectorField>& vectorFields,
    HashPtrTable<sphericalTensorField>& sphericalTensorFields,
    HashPtrTable<symmTensorField>& symmTensorFields,
    HashPtrTable<tensorField>& tensorFields
)
{
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().keyword() != "type" && iter().keyword() != "value")
        {
            if
            (
                iter().isStream()
             && iter().stream().size()
            )
            {
                ITstream& is = iter().stream();

                // Read first token
                token firstToken(is);

                if
                (
                    firstToken.isWord()
                 && firstToken.wordToken() == "nonuniform"
                )
                {
                    token fieldToken(is);

                    if (!fieldToken.isCompound())
                    {
                        if
                        (
                            fieldToken.isLabel()
                         && fieldToken.labelToken() == 0
                        )
                        {
                            // Ignore nonuniform 0 entry
                        }
                        else
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    token following 'nonuniform' "
                                  "is not a compound"
                            << exit(FatalIOError);
                        }
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<scalar>>::typeName
                    )
                    {
                        scalarField* fPtr = new scalarField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<scalar>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != patchSize_)
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << patchSize_ << ')'
                                << "\n    on patch " << patchName_
                                << exit(FatalIOError);
                        }

                        scalarFields.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<vector>>::typeName
                    )
                    {
                        vectorField* fPtr = new vectorField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<vector>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != patchSize_)
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << patchSize_ << ')'
                                << "\n    on patch " << patchName_
                                << exit(FatalIOError);
                        }

                        vectorFields.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<sphericalTensor>>::typeName
                    )
                    {
                        sphericalTensorField* fPtr = new sphericalTensorField;
                        fPtr->transfer
                        (
                            dynamicCast
                            <
                                token::Compound<List<sphericalTensor>>
                            >
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != patchSize_)
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << patchSize_ << ')'
                                << "\n    on patch " << patchName_
                                << exit(FatalIOError);
                        }

                        sphericalTensorFields.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<symmTensor>>::typeName
                    )
                    {
                        symmTensorField* fPtr = new symmTensorField;
                        fPtr->transfer
                        (
                            dynamicCast
                            <
                                token::Compound<List<symmTensor>>
                            >
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != patchSize_)
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << patchSize_ << ')'
                                << "\n    on patch " << patchName_
                                << exit(FatalIOError);
                        }

                        symmTensorFields.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<tensor>>::typeName
                    )
                    {
                        tensorField* fPtr = new tensorField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<tensor>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != patchSize_)
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << patchSize_ << ')'
                                << "\n    on patch " << patchName_
                                << exit(FatalIOError);
                        }

                        tensorFields.insert(iter().keyword(), fPtr);
                    }
                    else
                    {
                        FatalIOErrorInFunction
                        (
                            dict
                        )   << "\n    compound " << fieldToken.compoundToken()
                            << " not supported"
                            << "\n    on patch " << patchName_
                            << exit(FatalIOError);
                    }
                }
                else if
                (
                    firstToken.isWord()
                 && firstToken.wordToken() == "uniform"
                )
                {
                    token fieldToken(is);

                    if (!fieldToken.isPunctuation())
                    {
                        scalarFields.insert
                        (
                            iter().keyword(),
                            new scalarField
                            (
                                patchSize_,
                                fieldToken.number()
                            )
                        );
                    }
                    else
                    {
                        // Read as scalarList.
                        is.putBack(fieldToken);

                        scalarList l(is);

                        if (l.size() == vector::nComponents)
                        {
                            vector vs(l[0], l[1], l[2]);

                            vectorFields.insert
                            (
                                iter().keyword(),
                                new vectorField(patchSize_, vs)
                            );
                        }
                        else if (l.size() == sphericalTensor::nComponents)
                        {
                            sphericalTensor vs(l[0]);

                            sphericalTensorFields_.insert
                            (
                                iter().keyword(),
                                new sphericalTensorField(patchSize_, vs)
                            );
                        }
                        else if (l.size() == symmTensor::nComponents)
                        {
                            symmTensor vs(l[0], l[1], l[2], l[3], l[4], l[5]);

                            symmTensorFields.insert
                            (
                                iter().keyword(),
                                new symmTensorField(patchSize_, vs)
                            );
                        }
                        else if (l.size() == tensor::nComponents)
                        {
                            tensor vs
                            (
                                l[0], l[1], l[2],
                                l[3], l[4], l[5],
                                l[6], l[7], l[8]
                            );

                            tensorFields.insert
                            (
                                iter().keyword(),
                                new tensorField(patchSize_, vs)
                            );
                        }
                        else
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "\n    unrecognised native type " << l
                                << "\n    on patch " << patchName_
                                << exit(FatalIOError);
                        }
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::genericPatchFieldBase::genericPatchFieldBase
(
    const word& patchName,
    const label patchSize
)
:
    patchName_(patchName),
    patchSize_(patchSize)
{}


Foam::genericPatchFieldBase::genericPatchFieldBase
(
    const word& patchName,
    const label patchSize,
    const dictionary& dict
)
:
    patchName_(patchName),
    patchSize_(patchSize),
    actualTypeName_(dict.lookup("type")),
    dict_(dict)
{
    // Read contents into typed fields
    read
    (
        dict_,
        scalarFields_,
        vectorFields_,
        sphericalTensorFields_,
        symmTensorFields_,
        tensorFields_
    );
}


Foam::genericPatchFieldBase::genericPatchFieldBase
(
    const word& patchName,
    const label patchSize,
    const genericPatchFieldBase& ptf,
    const fvPatchFieldMapper& mapper
)
:
    patchName_(patchName),
    patchSize_(patchSize)
{
    // Write input field as a dictionary
    OStringStream os;
    os << token::BEGIN_BLOCK << nl;
    ptf.write(os);
    os << token::END_BLOCK << nl;
    IStringStream is(os.str());
    is >> dict_;

    actualTypeName_ = word(dict_.lookup("type"));

    // Read data fields from original patchField
    HashPtrTable<scalarField> scalarFields;
    HashPtrTable<vectorField> vectorFields;
    HashPtrTable<sphericalTensorField> sphericalTensorFields;
    HashPtrTable<symmTensorField> symmTensorFields;
    HashPtrTable<tensorField> tensorFields;
    read
    (
        dict_,
        scalarFields,
        vectorFields,
        sphericalTensorFields,
        symmTensorFields,
        tensorFields
    );



    // Map using mapper
    forAllConstIter
    (
        HashPtrTable<scalarField>,
        scalarFields,
        iter
    )
    {
        scalarFields_.insert
        (
            iter.key(),
            new scalarField(*iter(), mapper)
        );
    }

    forAllConstIter
    (
        HashPtrTable<vectorField>,
        vectorFields,
        iter
    )
    {
        vectorFields_.insert
        (
            iter.key(),
            new vectorField(*iter(), mapper)
        );
    }

    forAllConstIter
    (
        HashPtrTable<sphericalTensorField>,
        sphericalTensorFields,
        iter
    )
    {
        sphericalTensorFields_.insert
        (
            iter.key(),
            new sphericalTensorField(*iter(), mapper)
        );
    }

    forAllConstIter
    (
        HashPtrTable<symmTensorField>,
        symmTensorFields,
        iter
    )
    {
        symmTensorFields_.insert
        (
            iter.key(),
            new symmTensorField(*iter(), mapper)
        );
    }

    forAllConstIter
    (
        HashPtrTable<tensorField>,
        tensorFields,
        iter
    )
    {
        tensorFields_.insert
        (
            iter.key(),
            new tensorField(*iter(), mapper)
        );
    }
}


Foam::genericPatchFieldBase::genericPatchFieldBase
(
    const genericPatchFieldBase& ptf
)
:
    patchName_(ptf.patchName_),
    patchSize_(ptf.patchSize_),
    actualTypeName_(ptf.actualTypeName_),
    dict_(ptf.dict_),
    scalarFields_(ptf.scalarFields_),
    vectorFields_(ptf.vectorFields_),
    sphericalTensorFields_(ptf.sphericalTensorFields_),
    symmTensorFields_(ptf.symmTensorFields_),
    tensorFields_(ptf.tensorFields_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::genericPatchFieldBase::autoMap
(
    const fvPatchFieldMapper& m
)
{
    forAllIter
    (
        HashPtrTable<scalarField>,
        scalarFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }

    forAllIter
    (
        HashPtrTable<vectorField>,
        vectorFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }

    forAllIter
    (
        HashPtrTable<sphericalTensorField>,
        sphericalTensorFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }

    forAllIter
    (
        HashPtrTable<symmTensorField>,
        symmTensorFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }

    forAllIter
    (
        HashPtrTable<tensorField>,
        tensorFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }
}


void Foam::genericPatchFieldBase::rmap
(
    const genericPatchFieldBase& dptf,
    const labelList& addr
)
{
    // We're looping over the remote field entries so we can add any missing
    // 'nonuniform 0' entries

    forAllConstIter
    (
        HashPtrTable<scalarField>,
        dptf.scalarFields_,
        dptfIter
    )
    {
        HashPtrTable<scalarField>::iterator iter =
            scalarFields_.find(dptfIter.key());

        if (iter != scalarFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
        else
        {
            scalarField* ptr = new scalarField(patchSize_);
            scalarFields_.insert(dptfIter.key(), ptr);
            ptr->rmap(*dptfIter(), addr);
        }
    }

    forAllConstIter
    (
        HashPtrTable<vectorField>,
        dptf.vectorFields_,
        dptfIter
    )
    {
        HashPtrTable<vectorField>::const_iterator iter =
            vectorFields_.find(dptfIter.key());

        if (iter != vectorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
        else
        {
            vectorField* ptr = new vectorField(patchSize_);
            vectorFields_.insert(dptfIter.key(), ptr);
            ptr->rmap(*dptfIter(), addr);
        }
    }

    forAllConstIter
    (
        HashPtrTable<sphericalTensorField>,
        dptf.sphericalTensorFields_,
        dptfIter
    )
    {
        HashPtrTable<sphericalTensorField>::const_iterator iter =
            sphericalTensorFields_.find(dptfIter.key());

        if (iter != sphericalTensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
        else
        {
            sphericalTensorField* ptr = new sphericalTensorField(patchSize_);
            sphericalTensorFields_.insert(dptfIter.key(), ptr);
            ptr->rmap(*dptfIter(), addr);
        }
    }

    forAllConstIter
    (
        HashPtrTable<symmTensorField>,
        dptf.symmTensorFields_,
        dptfIter
    )
    {
        HashPtrTable<symmTensorField>::const_iterator iter =
            symmTensorFields_.find(dptfIter.key());

        if (iter != symmTensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
        else
        {
            symmTensorField* ptr = new symmTensorField(patchSize_);
            symmTensorFields_.insert(dptfIter.key(), ptr);
            ptr->rmap(*dptfIter(), addr);
        }
    }

    forAllConstIter
    (
        HashPtrTable<tensorField>,
        dptf.tensorFields_,
        dptfIter
    )
    {
        HashPtrTable<tensorField>::const_iterator iter =
            tensorFields_.find(dptfIter.key());

        if (iter != tensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
        else
        {
            tensorField* ptr = new tensorField(patchSize_);
            tensorFields_.insert(dptfIter.key(), ptr);
            ptr->rmap(*dptfIter(), addr);
        }
    }
}


Foam::List<Foam::Pair<Foam::word>>
Foam::genericPatchFieldBase::entryTypes() const
{
    label n =
        scalarFields_.size()
      + vectorFields_.size()
      + sphericalTensorFields_.size()
      + symmTensorFields_.size()
      + tensorFields_.size();

    List<Pair<word>> entries(n);
    n = 0;
    {
        const word& typeName = token::Compound<List<scalar>>::typeName;
        const wordList names(scalarFields_.sortedToc());
        forAll(names, i)
        {
            entries[n++] = Pair<word>(names[i], typeName);
        }
    }
    {
        const word& typeName = token::Compound<List<vector>>::typeName;
        const wordList names(vectorFields_.sortedToc());
        forAll(names, i)
        {
            entries[n++] = Pair<word>(names[i], typeName);
        }
    }
    {
        const word& typeName = token::Compound<List<sphericalTensor>>::typeName;
        const wordList names(sphericalTensorFields_.sortedToc());
        forAll(names, i)
        {
            entries[n++] = Pair<word>(names[i], typeName);
        }
    }
    {
        const word& typeName = token::Compound<List<symmTensor>>::typeName;
        const wordList names(symmTensorFields_.sortedToc());
        forAll(names, i)
        {
            entries[n++] = Pair<word>(names[i], typeName);
        }
    }
    {
        const word& typeName = token::Compound<List<tensor>>::typeName;
        const wordList names(tensorFields_.sortedToc());
        forAll(names, i)
        {
            entries[n++] = Pair<word>(names[i], typeName);
        }
    }
    return entries;
}


void Foam::genericPatchFieldBase::addEntry
(
    const word& key,
    const word& type
)
{
    if (type == token::Compound<List<scalar>>::typeName)
    {
        scalarFields_.insert(key, new scalarField(0));
    }
    else if (type == token::Compound<List<vector>>::typeName)
    {
        vectorFields_.insert(key, new vectorField(0));
    }
    else if (type == token::Compound<List<sphericalTensor>>::typeName)
    {
        sphericalTensorFields_.insert(key, new sphericalTensorField(0));
    }
    else if (type == token::Compound<List<symmTensor>>::typeName)
    {
        symmTensorFields_.insert(key, new symmTensorField(0));
    }
    else if (type == token::Compound<List<tensor>>::typeName)
    {
        tensorFields_.insert(key, new tensorField(0));
    }
    else
    {
        FatalErrorInFunction << "problem" << exit(FatalError);
    }
}


void Foam::genericPatchFieldBase::write(Ostream& os) const
{
    os.writeKeyword("type") << actualTypeName_ << token::END_STATEMENT << nl;

    forAllConstIter(dictionary, dict_, iter)
    {
        if (iter().keyword() != "type" && iter().keyword() != "value")
        {
            if
            (
                iter().isStream()
             && iter().stream().size()
             && iter().stream()[0].isWord()
             && iter().stream()[0].wordToken() == "nonuniform"
            )
            {
                if (scalarFields_.found(iter().keyword()))
                {
                    scalarFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
                else if (vectorFields_.found(iter().keyword()))
                {
                    vectorFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
                else if (sphericalTensorFields_.found(iter().keyword()))
                {
                    sphericalTensorFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
                else if (symmTensorFields_.found(iter().keyword()))
                {
                    symmTensorFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
                else if (tensorFields_.found(iter().keyword()))
                {
                    tensorFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
            }
            else
            {
               iter().write(os);
            }
        }
    }
}


Foam::label Foam::genericPatchFieldBase::findEntry
(
    const wordPairList& lst,
    const word& key
)
{
    forAll(lst, i)
    {
        if (lst[i].first() == key)
        {
            return i;
        }
    }
    return -1;
}


void Foam::genericPatchFieldBase::mergeEntries
(
    const word& actualTypeName,
    const wordPairList& entryTypes,
    HashTable<wordPairList>& patchFieldToTypes
)
{
    HashTable<wordPairList>::iterator fndType =
        patchFieldToTypes.find(actualTypeName);

    if (fndType == patchFieldToTypes.end())
    {
        patchFieldToTypes.insert(actualTypeName, entryTypes);
    }
    else
    {
        // Merge
        forAll(entryTypes, i)
        {
            const word& entry = entryTypes[i].first();

            label index = findEntry(fndType(), entry);
            if (index != -1)
            {
                if
                (
                    entryTypes[i].second()
                 != fndType()[index].second()
                )
                {
                    FatalErrorInFunction << "Inconsistent types"
                        << " for patchField " << actualTypeName
                        << exit(FatalError);
                }
            }
            else
            {
                fndType().append(entryTypes[i]);
            }
        }
    }
}


// ************************************************************************* //
