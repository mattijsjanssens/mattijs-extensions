/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2019 OpenCFD Ltd.
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

#include "FieldM.H"
#include "FieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func(Field<ReturnType>& res, const UList<Type>& f)                        \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_F(ReturnType, res, =, ::Foam::Func, Type, f)            \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func(const UList<Type>& f)                              \
{                                                                              \
    auto tres = tmp<Field<ReturnType>>::New(f.size());                         \
    Func(tres.ref(), f);                                                       \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func(const tmp<Field<Type>>& tf)                        \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type>::New(tf);                           \
    Func(tres.ref(), tf());                                                    \
    tf.clear();                                                                \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc(Field<ReturnType>& res, const UList<Type>& f)                      \
{                                                                              \
    TFOR_ALL_F_OP_OP_F(ReturnType, res, =, Op, Type, f)                        \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op(const UList<Type>& f)                       \
{                                                                              \
    auto tres = tmp<Field<ReturnType>>::New(f.size());                         \
    OpFunc(tres.ref(), f);                                                     \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op(const tmp<Field<Type>>& tf)                 \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type>::New(tf);                           \
    OpFunc(tres.ref(), tf());                                                  \
    tf.clear();                                                                \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_F_F                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, f1, Type2, f2                 \
    )                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    auto tres = tmp<Field<ReturnType>>::New(f1.size());                        \
    Func(tres.ref(), f1, f2);                                                  \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type2>::New(tf2);                         \
    Func(tres.ref(), f1, tf2());                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type1>::New(tf1);                         \
    Func(tres.ref(), tf1(), f2);                                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpTmp<ReturnType, Type1, Type1, Type2>::New(tf1, tf2);   \
    Func(tres.ref(), tf1(), tf2());                                            \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const Type1& s1,                                                           \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_S_F                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, s1, Type2, f2                 \
    )                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const Type1& s1,                                                           \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    auto tres = tmp<Field<ReturnType>>::New(f2.size());                        \
    Func(tres.ref(), s1, f2);                                                  \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type2>::New(tf2);                         \
    Func(tres.ref(), s1, tf2());                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_FUNC_F_S                                                     \
    (                                                                          \
        ReturnType, res, =, ::Foam::Func, Type1, f1, Type2, s2                 \
    )                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    auto tres = tmp<Field<ReturnType>>::New(f1.size());                        \
    Func(tres.ref(), f1, s2);                                                  \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> Func                                                    \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type1>::New(tf1);                         \
    Func(tres.ref(), tf1(), s2);                                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_F_OP_F(ReturnType, res, =, Type1, f1, Op, Type2, f2)         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op                                             \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    auto tres = tmp<Field<ReturnType>>::New(f1.size());                        \
    OpFunc(tres.ref(), f1, f2);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op                                             \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type2>::New(tf2);                         \
    OpFunc(tres.ref(), f1, tf2());                                             \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op                                             \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type1>::New(tf1);                         \
    OpFunc(tres.ref(), tf1(), f2);                                             \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op                                             \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpTmp<ReturnType, Type1, Type1, Type2>::New(tf1, tf2);   \
    OpFunc(tres.ref(), tf1(), tf2());                                          \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const Type1& s1,                                                           \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_S_OP_F(ReturnType, res, =, Type1, s1, Op, Type2, f2)         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op                                             \
(                                                                              \
    const Type1& s1,                                                           \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    auto tres = tmp<Field<ReturnType>>::New(f2.size());                        \
    OpFunc(tres.ref(), s1, f2);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op                                             \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<Field<Type2>>& tf2                                               \
)                                                                              \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type2>::New(tf2);                         \
    OpFunc(tres.ref(), s1, tf2());                                             \
    tf2.clear();                                                               \
    return tres;                                                               \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    Field<ReturnType>& res,                                                    \
    const UList<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    TFOR_ALL_F_OP_F_OP_S(ReturnType, res, =, Type1, f1, Op, Type2, s2)         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op                                             \
(                                                                              \
    const UList<Type1>& f1,                                                    \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    auto tres = tmp<Field<ReturnType>>::New(f1.size());                        \
    OpFunc(tres.ref(), f1, s2);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<Field<ReturnType>> operator Op                                             \
(                                                                              \
    const tmp<Field<Type1>>& tf1,                                              \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    auto tres = reuseTmp<ReturnType, Type1>::New(tf1);                         \
    OpFunc(tres.ref(), tf1(), s2);                                             \
    tf1.clear();                                                               \
    return tres;                                                               \
}


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)


// ************************************************************************* //
