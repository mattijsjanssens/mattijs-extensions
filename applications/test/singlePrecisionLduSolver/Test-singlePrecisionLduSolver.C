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




int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    


//    volScalarField fx(pow(mesh.C().component(vector::X), 1));
//    fx.write();
//    volScalarField gradx4(fvc::grad(fx)().component(vector::X));
//    gradx4.write();
//
//    volVectorField curlC(fvc::curl(1.0*mesh.C()));
//    curlC.write();
//
//    surfaceScalarField xf(mesh.Cf().component(vector::X));
//    surfaceScalarField xf4(pow(xf, 4));
//
//    for (int i=1; i<xf4.size()-1; i++)
//    {
//        scalar gradx4a = (xf4[i] - xf4[i-1])/(xf[i] - xf[i-1]);
//        Info<< (gradx4a - gradx4[i])/gradx4a << endl;
//    }

    Info<< "End" << endl;
}


// ************************************************************************* //
