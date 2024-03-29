/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Description
    Specialisation of Field\<T\> for scalar.

\*---------------------------------------------------------------------------*/

#include "scalarField.H"
#include "unitConversion.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
FieldWrapper<scalar, scalar>::FieldWrapper(scalarField& f)
:
    tmp<Field<scalar>>(f),
    ref_(f)
{}
template<>
FieldWrapper<solveScalar, solveScalar>::FieldWrapper(solveScalarField& f)
:
    tmp<Field<solveScalar>>(f),
    ref_(f)
{}
template<>
FieldWrapper<solveScalar, scalar>::FieldWrapper(scalarField& f)
:
    tmp<Field<solveScalar>>(tmp<Field<solveScalar>>::New(f.size())),
    ref_(f)
{
    copy(f, this->ref());
}
template<>
FieldWrapper<scalar, solveScalar>::FieldWrapper(solveScalarField& f)
:
    tmp<Field<scalar>>(tmp<Field<scalar>>::New(f.size())),
    ref_(f)
{
    copy(f, this->ref());
}

template<>
ConstFieldWrapper<scalar, scalar>::ConstFieldWrapper(const scalarField& f)
:
    tmp<Field<scalar>>(f)
{}
template<>
const scalarField&
ConstFieldWrapper<scalar, scalar>::get(const scalarField& f, scalarField&)
{
    return f;
}
template<>
ConstFieldWrapper<solveScalar, solveScalar>::ConstFieldWrapper
(
    const solveScalarField& f
)
:
    tmp<Field<solveScalar>>(f)
{}
template<>
const solveScalarField& ConstFieldWrapper<solveScalar, solveScalar>::
get(const solveScalarField& f, solveScalarField&)
{
    return f;
}

// Conversion specialisations
template<>
ConstFieldWrapper<solveScalar, scalar>::ConstFieldWrapper(const scalarField& f)
:
    tmp<Field<solveScalar>>(tmp<Field<solveScalar>>::New(f.size()))
{
    copy(f, this->ref());
}
template<>
const solveScalarField&
ConstFieldWrapper<solveScalar, scalar>::
get(const scalarField& f, solveScalarField& store)
{
    return copy(f, store);
}
template<>
ConstFieldWrapper<scalar, solveScalar>::ConstFieldWrapper
(
    const solveScalarField& f
)
:
    tmp<Field<scalar>>(tmp<Field<scalar>>::New(f.size()))
{
    copy(f, this->ref());
}
template<>
const scalarField&
ConstFieldWrapper<scalar, solveScalar>::
get(const solveScalarField& f, scalarField& store)
{
    return copy(f, store);
}


template<>
float solveVSmall<float>::value()
{
    return floatScalarVSMALL;
}
template<>
double solveVSmall<double>::value()
{
    return doubleScalarVSMALL;
}


template<>
tmp<scalarField> scalarField::component(const direction) const
{
    return *this;
}

void component(scalarField& sf, const UList<scalar>& f, const direction)
{
    sf = f;
}

template<>
void scalarField::replace(const direction, const UList<scalar>& sf)
{
    *this = sf;
}

template<>
void scalarField::replace(const direction, const scalar& s)
{
    *this = s;
}


void stabilise(scalarField& res, const UList<scalar>& sf, const scalar s)
{
    TFOR_ALL_F_OP_FUNC_S_F
    (
        scalar, res, =, ::Foam::stabilise, scalar, s, scalar, sf
    )
}

tmp<scalarField> stabilise(const UList<scalar>& sf, const scalar s)
{
    auto tresult = tmp<scalarField>::New(sf.size());
    stabilise(tresult.ref(), sf, s);
    return tresult;
}

tmp<scalarField> stabilise(const tmp<scalarField>& tsf, const scalar s)
{
    tmp<scalarField> tresult = New(tsf);
    stabilise(tresult.ref(), tsf(), s);
    tsf.clear();
    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
scalar sumProd(const UList<scalar>& f1, const UList<scalar>& f2)
{
    scalar SumProd = 0.0;
    if (f1.size() && (f1.size() == f2.size()))
    {
        TFOR_ALL_S_OP_F_OP_F(scalar, SumProd, +=, scalar, f1, *, scalar, f2)
    }
    return SumProd;
}


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


#define BesselFunc(func)                                                       \
void func(scalarField& res, const int n, const UList<scalar>& sf)              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_S_F(scalar, res, =, ::Foam::func, int, n, scalar, sf)   \
}                                                                              \
                                                                               \
tmp<scalarField> func(const int n, const UList<scalar>& sf)                    \
{                                                                              \
    auto tresult = tmp<scalarField>::New(sf.size());                           \
    func(tresult.ref(), n, sf);                                                \
    return tresult;                                                            \
}                                                                              \
                                                                               \
tmp<scalarField> func(const int n, const tmp<scalarField>& tsf)                \
{                                                                              \
    tmp<scalarField> tresult = New(tsf);                                       \
    func(tresult.ref(), n, tsf());                                             \
    tsf.clear();                                                               \
    return tresult;                                                            \
}

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
