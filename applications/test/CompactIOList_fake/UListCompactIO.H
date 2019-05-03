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

Class
    Foam::UListCompactIO

Description
    Compact storage of objects of type \<T\> (themselves List<BaseType>.

    Useful for lists of small sublists e.g. faceList, cellList.

SourceFiles
    UListCompactIO.C

\*---------------------------------------------------------------------------*/

#ifndef UListCompactIO_H
#define UListCompactIO_H

#include "labelList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
template<class T, class BaseType> class UListCompactIO;

template<class T, class BaseType>
void writeEntry(Ostream& os, const UListCompactIO<T, BaseType>& l);

template<class T, class BaseType> Istream& operator>>
(
    Istream&,
    UListCompactIO<T, BaseType>&
);

template<class T, class BaseType> Ostream& operator<<
(
    Ostream&,
    const UListCompactIO<T, BaseType>&
);


/*---------------------------------------------------------------------------*\
                           Class UListCompactIO Declaration
\*---------------------------------------------------------------------------*/

template<class T, class BaseType>
class UListCompactIO
{
protected:

        labelList start_;

        List<BaseType> elems_;

        // List of size+pointer. Should be UList but UList does not allow null
        // construct
        List<Tuple2<label, BaseType*>> ptrs_;


    // Protected Member Functions

        static void convertToCompact
        (
            labelList& start,
            List<BaseType>& elems,
            const UList<T>& lst
        );

        static void convertFromCompact
        (
            List<T>& lst,
            const labelUList& start,
            const UList<BaseType>& elems
        );

        static void setPointers
        (
            List<Tuple2<label, BaseType*>>& ptrs,
            const labelUList& start,
            UList<BaseType>& elems
        );


public:

    // Constructors

        //- Construct null
        UListCompactIO();

        //- Construct from UList
        UListCompactIO(const UList<T>&);

        //- Construct from Istream
        UListCompactIO(Istream&);
 

    // Member operators

        inline label size() const
        {
            return ptrs_.size();
        }

        //- Return element of UList. Only const access allowed since otherwise
        //  resizing the T would be allowed.
        inline const T& operator[](const label i) const
        {
            return reinterpret_cast<const T&>(ptrs_[i]);
        }


    // Member operators

        void operator=(const UListCompactIO<T, BaseType>&);

        void operator=(const UList<T>&);


    // IOstream functions

        friend void writeEntry <T, BaseType>
        (
            Ostream& os,
            const UListCompactIO<T, BaseType>& l
        );


    // IOstream operators

        //- Read List from Istream, discarding contents of existing List.
        friend Istream& operator>> <T, BaseType>
        (
            Istream&,
            UListCompactIO<T, BaseType>&
        );

        // Write List to Ostream.
        friend Ostream& operator<< <T, BaseType>
        (
            Ostream&,
            const UListCompactIO<T, BaseType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UListCompactIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
