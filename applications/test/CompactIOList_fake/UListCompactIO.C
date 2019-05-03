/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "UListCompactIO.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::UListCompactIO<T, BaseType>::convertToCompact
(
    labelList& start,
    List<BaseType>& elems,
    const UList<T>& lst
)
{
    start.setSize(lst.size() + 1);

    start[0] = 0;
    for (label i = 1; i < start.size(); i++)
    {
        label prev = start[i-1];
        start[i] = prev + lst[i-1].size();

        if (start[i] < prev)
        {
            FatalErrorInFunction
                << "Overall number of elements " << start[i]
                << " of UListCompactIO of size "
                << lst.size() << " overflows the representation of a label"
                << endl << "Please recompile with a larger representation"
                << " for label" << exit(FatalError);
        }
    }

    elems.setSize(start[start.size() - 1]);

    label elemi = 0;
    forAll(lst, i)
    {
        const T& subList = lst[i];

        forAll(subList, j)
        {
            elems[elemi++] = subList[j];
        }
    }
}


template<class T, class BaseType>
void Foam::UListCompactIO<T, BaseType>::convertFromCompact
(
    List<T>& lst,
    const labelUList& start,
    const UList<BaseType>& elems
)
{
    lst.setSize(start.size() - 1);

    forAll(lst, i)
    {
        T& subList = lst[i];

        label index = start[i];
        subList.setSize(start[i+1] - index);

        forAll(subList, j)
        {
            subList[j] = elems[index++];
        }
    }
}


template<class T, class BaseType>
void Foam::UListCompactIO<T, BaseType>::setPointers
(
    List<Tuple2<label, BaseType*>>& ptrs,
    const labelUList& start,
    UList<BaseType>& elems
)
{
    ptrs.setSize(start.size()-1);
    forAll(ptrs, i)
    {
        const label sz = start[i+1]-start[i];
        ptrs[i].first() = sz;
        ptrs[i].second() = &(elems[start[i]]);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::UListCompactIO<T, BaseType>::UListCompactIO()
{}


template<class T, class BaseType>
Foam::UListCompactIO<T, BaseType>::UListCompactIO(const UList<T>& l)
:
    ptrs_(l.size())
{
    convertToCompact(start_, elems_, l);
    setPointers(ptrs_, start_, elems_);
}


template<class T, class BaseType>
Foam::UListCompactIO<T, BaseType>::UListCompactIO(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::UListCompactIO<T, BaseType>::operator=
(
    const UListCompactIO<T, BaseType>& rhs
)
{
    start_ = rhs.start_;
    elems_ = rhs.elems_;
    setPointers(ptrs_, start_, elems_);
}


template<class T, class BaseType>
void Foam::UListCompactIO<T, BaseType>::operator=(const UList<T>& rhs)
{
    convertToCompact(start_, elems_, rhs);
    setPointers(ptrs_, start_, elems_);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::writeEntry(Ostream& os, const UListCompactIO<T, BaseType>& l)
{
    writeEntry(os, l.start_);
    writeEntry(os, l.elems_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::UListCompactIO<T, BaseType>& L
)
{
    is >> L.start_ >> L.elems_;
    UListCompactIO<T, BaseType>::setPointers(L.ptrs_, L.start_, L.elems_);

    return is;
}


template<class T, class BaseType>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::UListCompactIO<T, BaseType>& L
)
{
    os << L.start_ << L.elems_;

    return os;
}


// ************************************************************************* //
