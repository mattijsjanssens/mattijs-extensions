/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    test

Description
    Finite volume method test code.

\*---------------------------------------------------------------------------*/

#include "scalarField.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//const scalarField& get(const solveScalarField& f, scalarField& store)
//{
//    #ifdef WM_DP
//    return f;
//    #else
//    store.setSize(f.size());
//    forAll(store, i)
//    {
//        store[i] = f[i];
//    }
//    return store;
//    #endif
//}
//tmp<scalarField> get(const solveScalarField& f)
//{
//    #ifdef WM_DP
//    return f;
//    #else
//    tmp<scalarField> tstore(tmp<scalarField>::New(f.size());
//    scalarField& store = tstore.ref();
//    forAll(store, i)
//    {
//        store[i] = f[i];
//    }
//    return tstore;
//    #endif
//}


//
//
//class solveScalarField
//:
//    public Field<solveScalar>
//{
//public:
//
//    // Constructors
//
//        //- Construct from scalarField
//        solveScalarWrapper(const scalarField& fld)
//        :
//            Field<solveScalar>(fld.size())
//        {
//            forAll(fld, i)
//            {
//                operator[](i] = fld[i];
//            }
//        }
//    {}
//}
//
//
//
//class solveScalarWrapper
//{
//    // Private data
//
//        //- Reference to scalarField
//        tmp<scalarField> scalars_;
//
//        //- Reference to solveScalarField
//        tmp<solveScalarField> solveScalars_;
//
//
//    // Private Member Functions
//
//        //- No copy construct
//        solveScalarWrapper(const solveScalarWrapper&) = delete;
//
//        //- No copy assignment
//        void operator=(const solveScalarWrapper&) = delete;
//
//
//public:
//
//    // Constructors
//
//        //- Construct from components
//        solveScalarWrapper(scalarField& data, solveScalarField& solveData)
//        :
//            data_(data),
//            solveData_(solveData)
//        {}
//
//
////    //- Destructor
////    ~solveScalarWrapper();
//
//
//    // Member Functions
//
//        // Access
//
//            solveScalarField& ref()
//            {
//                #ifdef WM_DP
//                    return data_;
//                #else
//                    if (
//                #endif
//            }
//
//
//        // Check
//
//        // Edit
//
//        // Write
//
//
//    // Member Operators
//
//        void operator=(const solveScalarWrapper&);
//
//
//    // Friend Functions
//
//    // Friend Operators
//
//    // IOstream Operators
//
//        friend Istream& operator>>(Istream&, solveScalarWrapper&);
//        friend Ostream& operator<<(Ostream&, const solveScalarWrapper&);
//};
//
//
//// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
//} // End namespace Foam
//
//// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
//#include "solveScalarWrapperI.H"
//
//// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
//#endif


typedef Field<scalar> scalarField;
typedef Field<solveScalar> solveScalarField;

template<class Type, class InputType>
class FieldWrapper2
:
    public tmp<Field<Type>>
{
    // Private data

        //- Reference to underlying field
        Field<InputType>& ref_;


public:

    // Constructors

        //- Construct from Field<InputType>
        FieldWrapper2(Field<InputType>& f) = delete;

    //- Destructor
    ~FieldWrapper2()
    {
        if (this->isTmp())
        {
            const Field<Type>& store = this->operator()();
            ref_.setSize(store.size());
            forAll(ref_, i)
            {
                ref_[i] = store[i];
            }
        }
    }
};

template<>
FieldWrapper2<scalar, scalar>::FieldWrapper2(scalarField& f)
:
    tmp<Field<scalar>>(f),
    ref_(f)
{}

template<>
FieldWrapper2<solveScalar, scalar>::FieldWrapper2(scalarField& f)
:
    tmp<Field<solveScalar>>(tmp<Field<solveScalar>>::New(f.size())),
    ref_(f)
{
    Field<solveScalar>& store = this->ref();
    forAll(store, i)
    {
        store[i] = f[i];
    }
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    //#include "createMesh.H"

    scalarField input(3);
    FieldWrapper2<solveScalar, scalar> toutput(input);
    solveScalarField& output = toutput.constCast();

    Info<< "End" << endl;
}


// ************************************************************************* //
