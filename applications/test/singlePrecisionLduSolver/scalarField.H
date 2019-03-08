/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Typedef
    Foam::scalarField

Description
    Specialisation of Field\<T\> for scalar.

SourceFiles
    scalarField.C

\*---------------------------------------------------------------------------*/

#ifndef scalarField_H
#define scalarField_H

#include "Field.H"
#include "scalar.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<scalar> scalarField;
typedef Field<solveScalar> solveScalarField;


//- A Field wrapper with possible data conversion
template<class Type, class InputType>
class FieldWrapper
:
    public tmp<Field<Type>>
{
    // Private Data

        //- Reference to underlying field
        Field<InputType>& ref_;


public:

    // Constructors

        //- Construct from Field<InputType>, copying on input as required
        FieldWrapper(Field<InputType>& input)
        :
            tmp<Field<Type>>(),
            ref_(input)
        {
            if (std::is_same<Type, InputType>::value)
            {
                // Future
                #if 0
                this->cref(input);
                #else
                tmp<Field<Type>> other(reinterpret_cast<Field<Type>&>(input));
                this->swap(other);
                #endif
            }
            else
            {
                this->reset(new Field<Type>(input.size()));
                std::copy(input.cbegin(), input.cend(), this->ref().begin());
            }
        }


    //- Destructor, copying on destroy
    ~FieldWrapper()
    {
        if (this->isTmp())
        {
            const Field<Type>& store = this->operator()();
            ref_.setSize(store.size());
            std::copy(store.cbegin(), store.cend(), ref_.begin());
          
        }
    }
};


//- A Field wrapper with possible data conversion
template<class Type, class InputType>
class ConstFieldWrapper
:
    public tmp<Field<Type>>
{
public:

    // Constructors

        //- Construct from InputType
        ConstFieldWrapper(const Field<InputType>& input)
        :
            tmp<Field<Type>>()
        {
            if (std::is_same<Type, InputType>::value)
            {
                // Future
                #if 0
                this->cref(input);
                #else
                tmp<Field<Type>> other
                (
                    reinterpret_cast<const Field<Type>&>(input)
                );
                this->swap(other);
                #endif
            }
            else
            {
                this->reset(new Field<Type>(input.size()));
                std::copy(input.cbegin(), input.cend(), this->ref().begin());
            }
        }


    // Member Functions

       //- Return the field

        static const Field<Type>& get
        (
            const Field<InputType>& input,
            Field<Type>& dst
        )
        {
            if (std::is_same<Type, InputType>::value)
            {
                return reinterpret_cast<const Field<Type>&>(input);
            }
            else
            {
                dst.setSize(input.size());
                std::copy(input.cbegin(), input.cend(), dst.begin());
                return dst;
            }
        }
};



template<class Type>
class solveVSmall
{
public:

    static Type value() = delete;
};
template<>
float solveVSmall<float>::value();
template<>
double solveVSmall<double>::value();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<scalarField> scalarField::component(const direction) const;

void component
(
    scalarField& sf,
    const UList<scalar>& f,
    const direction
);

template<>
void scalarField::replace(const direction, const UList<scalar>& sf);

template<>
void scalarField::replace(const direction, const scalar& s);


void stabilise(scalarField& Res, const UList<scalar>& sf, const scalar s);
tmp<scalarField> stabilise(const UList<scalar>&, const scalar s);
tmp<scalarField> stabilise(const tmp<scalarField>&, const scalar s);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
float sumProd(const UList<float>& f1, const UList<float>& f2);

template<>
double sumProd(const UList<double>& f1, const UList<double>& f2);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, divide)

BINARY_FUNCTION(scalar, scalar, scalar, pow)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, pow)

BINARY_FUNCTION(scalar, scalar, scalar, atan2)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, atan2)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, scalar, pow3)
UNARY_FUNCTION(scalar, scalar, pow4)
UNARY_FUNCTION(scalar, scalar, pow5)
UNARY_FUNCTION(scalar, scalar, pow6)
UNARY_FUNCTION(scalar, scalar, pow025)
UNARY_FUNCTION(scalar, scalar, sqrt)
UNARY_FUNCTION(scalar, scalar, cbrt)
UNARY_FUNCTION(scalar, scalar, sign)
UNARY_FUNCTION(scalar, scalar, pos)
UNARY_FUNCTION(scalar, scalar, pos0)
UNARY_FUNCTION(scalar, scalar, neg)
UNARY_FUNCTION(scalar, scalar, neg0)
UNARY_FUNCTION(scalar, scalar, posPart)
UNARY_FUNCTION(scalar, scalar, negPart)
UNARY_FUNCTION(scalar, scalar, exp)
UNARY_FUNCTION(scalar, scalar, log)
UNARY_FUNCTION(scalar, scalar, log10)
UNARY_FUNCTION(scalar, scalar, sin)
UNARY_FUNCTION(scalar, scalar, cos)
UNARY_FUNCTION(scalar, scalar, tan)
UNARY_FUNCTION(scalar, scalar, asin)
UNARY_FUNCTION(scalar, scalar, acos)
UNARY_FUNCTION(scalar, scalar, atan)
UNARY_FUNCTION(scalar, scalar, sinh)
UNARY_FUNCTION(scalar, scalar, cosh)
UNARY_FUNCTION(scalar, scalar, tanh)
UNARY_FUNCTION(scalar, scalar, asinh)
UNARY_FUNCTION(scalar, scalar, acosh)
UNARY_FUNCTION(scalar, scalar, atanh)
UNARY_FUNCTION(scalar, scalar, erf)
UNARY_FUNCTION(scalar, scalar, erfc)
UNARY_FUNCTION(scalar, scalar, lgamma)
UNARY_FUNCTION(scalar, scalar, j0)
UNARY_FUNCTION(scalar, scalar, j1)
UNARY_FUNCTION(scalar, scalar, y0)
UNARY_FUNCTION(scalar, scalar, y1)

UNARY_FUNCTION(scalar, scalar, degToRad)
UNARY_FUNCTION(scalar, scalar, radToDeg)
UNARY_FUNCTION(scalar, scalar, atmToPa)
UNARY_FUNCTION(scalar, scalar, paToAtm)

#define BesselFunc(func)                                            \
void func(scalarField& Res, const int n, const UList<scalar>& sf);  \
tmp<scalarField> func(const int n, const UList<scalar>&);           \
tmp<scalarField> func(const int n, const tmp<scalarField>&);

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
