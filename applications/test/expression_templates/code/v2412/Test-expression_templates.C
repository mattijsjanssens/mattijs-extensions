/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 M. Janssens
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
    Test-expression_templates.C

Description
    Experiment with expression templates.
    See https://en.wikipedia.org/wiki/Expression_templates

\*---------------------------------------------------------------------------*/


#include <cassert>
#include "argList.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "Time.H"
#include "OBJstream.H"
#include "DynamicField.H"
#include "cyclicAMIPolyPatch.H"
#include "columnFvMesh.H"
#include "uindirectPrimitivePatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <typename E>
class VecExpression
{
public:
    static constexpr bool is_leaf = false;

    double operator[](size_t i) const
    {
        // Delegation to the actual expression type. This avoids dynamic
        // polymorphism (a.k.a. virtual functions in C++)
        return static_cast<E const&>(*this)[i];
    }
    size_t size() const { return static_cast<E const&>(*this).size(); }
};


class Vec
:
    public VecExpression<Vec>
{
    std::array<double, 3> elems;

public:
    static constexpr bool is_leaf = true;

    decltype(auto) operator[](size_t i) const
    {
        return elems[i];
    }
//    decltype(auto) &operator[](size_t i)
//    {
//        return elems[i];
//    }
    size_t size() const
    {
        return elems.size();
    }

    // construct Vec using initializer list 
    Vec(std::initializer_list<double> init)
    {
        std::copy(init.begin(), init.end(), elems.begin());
    }

    // A Vec can be constructed from any VecExpression, forcing its evaluation.
    template <typename E>
    Vec(VecExpression<E> const& expr)
    {
        for (size_t i = 0; i != expr.size(); ++i)
        {
            elems[i] = expr[i];
        }
    }
};


// Sum
// ~~~
// wrapper class + operator

template <typename E1, typename E2>
class VecSum
:
    public VecExpression<VecSum<E1, E2> >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type _u;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type _v;

public:
    static constexpr bool is_leaf = false;

    VecSum(E1 const& u, E2 const& v)
    :
     _u(u), _v(v)
    {
        assert(u.size() == v.size());
    }
    decltype(auto) operator[](size_t i) const
    {
        return _u[i] + _v[i];
    }
    size_t size() const { return _v.size(); }
};

template <typename E1, typename E2>
VecSum<E1, E2>
operator+(VecExpression<E1> const& u, VecExpression<E2> const& v)
{
    return VecSum<E1, E2>
    (
        *static_cast<const E1*>(&u),
        *static_cast<const E2*>(&v)
    );
}


// Diff
// ~~~~
// wrapper class + operator

template <typename E1, typename E2>
class VecDiff
:
    public VecExpression<VecDiff<E1, E2> >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type _u;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type _v;

public:
    static constexpr bool is_leaf = false;

    VecDiff(E1 const& u, E2 const& v)
    :
     _u(u), _v(v)
    {
        assert(u.size() == v.size());
    }
    decltype(auto) operator[](size_t i) const
    {
        return _u[i] - _v[i];
    }
    size_t size() const { return _v.size(); }
};

template <typename E1, typename E2>
VecDiff<E1, E2>
operator-(VecExpression<E1> const& u, VecExpression<E2> const& v)
{
    return VecDiff<E1, E2>
    (
        *static_cast<const E1*>(&u),
        *static_cast<const E2*>(&v)
    );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Try fields

template <typename E>
class FieldExpression
{
public:
    static constexpr bool is_leaf = false;

    double operator[](size_t i) const
    {
        // Delegation to the actual expression type. This avoids dynamic
        // polymorphism (a.k.a. virtual functions in C++)
        return static_cast<E const&>(*this)[i];
    }
    size_t size() const { return static_cast<E const&>(*this).size(); }
};


class FieldWrap
:
    public FieldExpression<FieldWrap>
{
    scalarField elems;

public:
    static constexpr bool is_leaf = true;

//    operator const scalarField&() const { return elems; }

    decltype(auto) operator[](size_t i) const
    {
        return elems[i];
    }
//    decltype(auto) &operator[](size_t i)
//    {
//        return elems[i];
//    }
    size_t size() const
    {
        return elems.size();
    }

    // construct FieldWrap using initializer list 
    FieldWrap(std::initializer_list<double> init)
    :
        elems(init)
    {}

    // A FieldWrap can be constructed from any FieldExpression, forcing its evaluation.
    template <typename E>
    FieldWrap(FieldExpression<E> const& expr)
    {
        for (size_t i = 0; i != expr.size(); ++i)
        {
            elems[i] = expr[i];
        }
    }
};


// Sum
// ~~~
// wrapper class + operator

template <typename E1, typename E2>
class FieldWrapSum
:
    public FieldExpression<FieldWrapSum<E1, E2> >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type _u;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type _v;

public:
    static constexpr bool is_leaf = false;

    FieldWrapSum(E1 const& u, E2 const& v)
    :
     _u(u), _v(v)
    {
        assert(u.size() == v.size());
    }
    decltype(auto) operator[](size_t i) const
    {
        return _u[i] + _v[i];
    }
    size_t size() const { return _v.size(); }
};

template <typename E1, typename E2>
FieldWrapSum<E1, E2>
operator+(FieldExpression<E1> const& u, FieldExpression<E2> const& v)
{
    return FieldWrapSum<E1, E2>
    (
        *static_cast<const E1*>(&u),
        *static_cast<const E2*>(&v)
    );
}


// Diff
// ~~~~
// wrapper class + operator

template <typename E1, typename E2>
class FieldWrapDiff
:
    public FieldExpression<FieldWrapDiff<E1, E2> >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type _u;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type _v;

public:
    static constexpr bool is_leaf = false;

    FieldWrapDiff(E1 const& u, E2 const& v)
    :
     _u(u), _v(v)
    {
        assert(u.size() == v.size());
    }
    decltype(auto) operator[](size_t i) const
    {
        return _u[i] - _v[i];
    }
    size_t size() const { return _v.size(); }
};

template <typename E1, typename E2>
FieldWrapDiff<E1, E2>
operator-(FieldExpression<E1> const& u, FieldExpression<E2> const& v)
{
    return FieldWrapDiff<E1, E2>
    (
        *static_cast<const E1*>(&u),
        *static_cast<const E2*>(&v)
    );
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



// Main
// ~~~~

int main() {
    FieldWrap v0({23.4,  12.5,  144.56});
    FieldWrap v1({67.12, 34.8,  90.34});
    FieldWrap v2({34.90, 111.9, 45.12});
    
//    // Following assignment will call the ctor of Vec which accept type of 
//    // `VecExpression<E> const&`. Then expand the loop body to 
//    // a.elems[i] + b.elems[i] + c.elems[i]
//    Vec sum_of_vec_type = v0 + v1 + v2; 
//
//    for (size_t i = 0; i < sum_of_vec_type.size(); ++i)
//        std::cout << sum_of_vec_type[i] << std::endl;

    // To avoid creating any extra storage, other than v0, v1, v2
    // one can do the following (Tested with C++11 on GCC 5.3.0)
    auto sum = v0 + v1 - v2;

    for (size_t i = 0; i < sum.size(); ++i)
        std::cout << sum[i] << std::endl;
    // Observe that in this case typeid(sum) will be VecSum<VecSum<Vec, Vec>, Vec>
    // and this chaining of operations can go on.

    scalarField result(sum.size());
    for (size_t i = 0; i < sum.size(); ++i)
    {
        result[i] = sum[i];
    }
    Pout<< "result:" << result << endl;


    return 0;
}


// ************************************************************************* //
