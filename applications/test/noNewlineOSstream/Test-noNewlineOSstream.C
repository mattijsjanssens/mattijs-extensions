/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Application
    Test-noNewlineOSstream

Description
    Multi-pass printing:

    mpOSstream os;
    os << dict;
    Pout<< os.str() << endl; // but prints with ""


    Or:

    streamDecorator os(Pout);
    os << dict;

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "noNewlineOSstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Container>
label countChars(const Container& ll)
{
    label sz = 0;
    forAll(ll, i)
    {
        sz += countChars(ll[i]) + 1;
    }
    return sz;
}
// Specialisation for labelList
template<>
label countChars(const labelList& ll)
{
    return ll.size()+1;
}


template<class Container>
void printWithLines(const Container& l)
{
    const label nChars = countChars(l);
    if (nChars > 10)
    {
        // Write size and start delimiter
        Pout<< l.size() << nl << indent << token::BEGIN_LIST;

        // Write contents
        if (l.size())
        {
            Pout<< incrIndent;
            forAll(l, i)
            {
                Pout<< nl << indent;
                printWithLines(l[i]);
            }
            Pout<< decrIndent << nl << indent << token::END_LIST;
        }
        else
        {
            // Write end delimiter
            Pout<< token::END_LIST;
        }
    }
    else
    {
        // Write size and start delimiter
        Pout<< l.size() << token::BEGIN_LIST;

        // Write contents
        forAll(l, i)
        {
            if (i > 0) Pout<< token::SPACE;
            Pout<< l[i];
        }

        // Write end delimiter
        Pout<< token::END_LIST;
    }
}
// Specialisation for label
template<>
void printWithLines(const label& l)
{
    Pout << l;
}



void indent(const label indent, const string& s, label& index)
{
    if (s[index] == '(')
    {
        
}



int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    labelList a(identity(11));

//     {
//         Pout<< "a:" << (noNewlineOSstream(Pout) << a);
//     }

    Pout<< "here" << endl;

//const noNewlineOSstream Pout2(Pout);
//Pout2 << "a:" << a << endl;


{
/*

((1 3 4)(8 9 7 8))

Count chars if printed on single line:
    18
If above limit decide to break:
    (
        (1 3 4)
        (8 9 7 8 5)
    )
or
    (
        (1 3 4)
        (
            8
            9
            7
            8
            5
        )
    )
or even
    (
        (1 3 4)
        (
            8 9
            7 8
            5
        )
    )
*/
    labelListList ll
    (
        {
            {1, 3, 4},
            {8, 9, 7, 8, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1}
        }
    );
    printWithLines(ll);
}

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
