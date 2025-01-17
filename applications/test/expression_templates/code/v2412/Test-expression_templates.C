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

    FieldWrap(const FieldWrap&) = default;

    // returns the underlying data
    const scalarField& data() const
    { 
        return elems; 
    }

    scalar operator[](size_t i) const
    {
        return elems[i];
    }
    scalar& operator[](size_t i)
    {
        return elems[i];
    }
    size_t size() const
    {
        return elems.size();
    }

    // construct FieldWrap using initializer list 
    FieldWrap(std::initializer_list<double> init)
    :
        elems(init)
    {}

    // A FieldWrap can be constructed from any FieldExpression, forcing its
    // evaluation.
    template <typename E>
    FieldWrap(FieldExpression<E> const& expr)
    {
        for (size_t i = 0; i != expr.size(); ++i)
        {
            elems[i] = expr[i];
        }
    }
};
class FieldRefWrap
:
    public FieldExpression<FieldRefWrap>
{
    const scalarField& elems;

public:
    // ! Store as copy (since holds reference)
    static constexpr bool is_leaf = false;

    // returns the underlying data
    const scalarField& data() const
    { 
        return elems; 
    }

    scalar operator[](size_t i) const
    {
        return elems[i];
    }
    size_t size() const
    {
        return elems.size();
    }

    // construct FieldRefWrap using initializer list 
    FieldRefWrap(const scalarField& init)
    :
        elems(init)
    {}
};


// Sum
// ~~~
// wrapper class + operator

template <typename E1, typename E2>
class FieldSum
:
    public FieldExpression<FieldSum<E1, E2> >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type _u;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type _v;

public:
    static constexpr bool is_leaf = false;

    FieldSum(E1 const& u, E2 const& v)
    :
     _u(u), _v(v)
    {
        assert(u.size() == v.size());
    }
    decltype(auto) operator[](size_t i) const
    {
        return _u[i] + _v[i];
    }
    size_t size() const { return _u.size(); }
};
template <typename E1, typename E2>
FieldSum<E1, E2>
operator+(FieldExpression<E1> const& u, FieldExpression<E2> const& v)
{
    return FieldSum<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}


// Diff
// ~~~~
// wrapper class + operator

template <typename E1, typename E2>
class FieldDiff
:
    public FieldExpression<FieldDiff<E1, E2> >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type _u;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type _v;

public:
    static constexpr bool is_leaf = false;

    FieldDiff(E1 const& u, E2 const& v)
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
FieldDiff<E1, E2>
operator-(FieldExpression<E1> const& u, FieldExpression<E2> const& v)
{
    return FieldDiff<E1, E2>
    (
        *static_cast<const E1*>(&u),
        *static_cast<const E2*>(&v)
    );
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template <typename E, typename WrapType>
class VolFieldExpression
{
public:
    static constexpr bool is_leaf = false;

    typedef WrapType wrapType;

    double operator[](size_t i) const
    {
        // Delegation to the actual expression type. This avoids dynamic
        // polymorphism (a.k.a. virtual functions in C++)
        return static_cast<E const&>(*this)[i];
    }
    size_t size() const { return static_cast<E const&>(*this).size(); }
};


class volScalarFieldRefWrap
:
    public VolFieldExpression<volScalarFieldRefWrap, FieldRefWrap>
{
    const volScalarField& elems;

public:

    static constexpr bool is_leaf = true;

    //- Type to return for patchField
    typedef FieldRefWrap wrapType;

    // returns the underlying data
    const volScalarField& data() const
    { 
        return elems; 
    }

    scalar operator[](size_t i) const
    {
        return elems[i];
    }
    size_t size() const
    {
        return elems.size();
    }

    wrapType patchField(const label i) const
    {
        return elems.boundaryField()[i];
    }

    // construct volScalarFieldRefWrap from components
    volScalarFieldRefWrap(const volScalarField& init)
    :
        elems(init)
    {}
};


template <typename E1, typename E2>
class volScalarFieldSum
:
    public VolFieldExpression
    <
        volScalarFieldSum<E1, E2>,
        FieldSum
        <
            typename E1::wrapType,
            typename E2::wrapType
        >
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type _u;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type _v;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for patchField
    typedef FieldSum
    <
        typename E1::wrapType,
        typename E2::wrapType
    > wrapType;

    volScalarFieldSum(E1 const& u, E2 const& v)
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

    wrapType patchField(const label i) const
    {
        return wrapType(_u.patchField(i), _v.patchField(i));
    }
};
template <typename E1, typename E2>
volScalarFieldSum<E1, E2>
operator+
(
    VolFieldExpression<E1, typename E1::wrapType> const& u,
    VolFieldExpression<E2, typename E2::wrapType> const& v
)
{
    return volScalarFieldSum<E1, E2>
    (
        *static_cast<const E1*>(&u),
        *static_cast<const E2*>(&v)
    );
}
template <typename E1, typename E2>
class volScalarFieldDiff
:
    public VolFieldExpression
    <
        volScalarFieldDiff<E1, E2>,
        FieldDiff
        <
            typename E1::wrapType,
            typename E2::wrapType
        >
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type _u;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type _v;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for patchField
    typedef FieldDiff
    <
        typename E1::wrapType,
        typename E2::wrapType
    > wrapType;

    volScalarFieldDiff(E1 const& u, E2 const& v)
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

    wrapType patchField(const label i) const
    {
        return wrapType(_u.patchField(i), _v.patchField(i));
    }
};
template <typename E1, typename E2>
volScalarFieldDiff<E1, E2>
operator-
(
    VolFieldExpression<E1, typename E1::wrapType> const& u,
    VolFieldExpression<E2, typename E2::wrapType> const& v
)
{
    return volScalarFieldDiff<E1, E2>
    (
        *static_cast<const E1*>(&u),
        *static_cast<const E2*>(&v)
    );
}


// Main
// ~~~~

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

if (false)
{
    FieldWrap v0({23.4,  12.5,  144.56});
    FieldWrap v1({67.12, 34.8,  90.34});
    FieldWrap v2({34.90, 111.9, 45.12});
    
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
}

if (false)
{
    const scalarField v0Data({23.4,  12.5,  144.56});
    const scalarField v1Data({67.12, 34.8,  90.34});
    const scalarField v2Data({34.90, 111.9, 45.12});
    
    const FieldRefWrap v0(v0Data);
    const FieldRefWrap v1(v1Data);
    const FieldRefWrap v2(v2Data);
    
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
    Pout<< "result2:" << result << endl;
}
{
    const volScalarField fld0
    (
        IOobject
        (
            "dummy0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );
    const volScalarField fld1
    (
        IOobject
        (
            "dummy1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );
    const volScalarField fld2
    (
        IOobject
        (
            "dummy2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    const volScalarFieldRefWrap wfld0(fld0);
    const volScalarFieldRefWrap wfld1(fld1);
    const volScalarFieldRefWrap wfld2(fld2);

    
    // To avoid creating any extra storage, other than v0, v1, v2
    // one can do the following (Tested with C++11 on GCC 5.3.0)
    auto sum = wfld0 + wfld1 - wfld2;

    for (size_t i = 0; i < sum.size(); ++i)
    {
        std::cout << sum[i] << std::endl;
    }

    const label n = fld0.boundaryField().size();

    for (label i = 0; i < n; ++i)
    {
        Pout<< "Patch:" << i << endl;

        // 1. Repeat expression : works
        //auto patchSum =
        //    FieldRefWrap(fld0.boundaryField()[i])
        //  + FieldRefWrap(fld1.boundaryField()[i]);


        // 2. Use original class member function
        auto patchSum = sum.patchField(i);

        Pout<< "Patch size:" << patchSum.size() << endl;

        for (size_t i = 0; i < patchSum.size(); ++i)
        {
            std::cout << patchSum[i] << std::endl;
        }
    }
}
    return 0;
}


// ************************************************************************* //
