/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Mattijs Janssens
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

Description:
    See
    https://stackoverflow.com/questions/7868936/read-file-line-by-line-using-ifstream-in-c

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "argList.H"
#include "Time.H"
#include "OFstream.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// void write(Ostream& os, const UList<label>& L)
// {
//     if (os.format() == IOstream::ASCII)
//     {
//         // Write size and start delimiter
//         os << nl << L.size() << nl << token::BEGIN_LIST;
// 
//     }


// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    //#include "createTime.H"
    //#include "createMesh.H"

    scalarField fld(10000000, 0.0);
    forAll(fld, i)
    {
        fld[i] = 1.0/(i+1);
    }

    // Write scalarField in ascii
    cpuTime timer;
    {
        OFstream os("scalars_foam.txt");
        os << fld << endl;
    }
    Pout<< "Written OpenFOAM " << fld.size() << " in "
        << timer.cpuTimeIncrement()
        << endl;
        
    {
        FILE* fp = fopen("scalars_c.txt", "w");
        forAll(fld, i)
        {
            fprintf(fp, "%g\n", fld[i]);
        }
        fclose(fp);
    }
    Pout<< "Written C " << fld.size() << " in " << timer.cpuTimeIncrement()
        << endl;


    // Read
    {
        std::ifstream file("scalars_c.txt");
        if (file.is_open())
        {
            label i = 0;
            std::string line;
            while (getline(file, line))
            {
                sscanf(line.c_str(), "%lg", &fld[i]);
                i++;
            }
            file.close();
        }
    } 
    Pout<< "Read C++ " << fld.size() << " in " << timer.cpuTimeIncrement()
        << endl;

    {
        IFstream file("scalars_c.txt");
        if (file.good())
        {
            label i = 0;
            while (file.good())
            {
                file >> fld[i];
                i++;
            }
        }
    } 
    Pout<< "Read OpenFOAM " << fld.size() << " in " << timer.cpuTimeIncrement()
        << endl;

    {
        FILE* fp = fopen("scalars_c.txt", "r");
        if (fp == NULL)
            exit(1);

        char* line = NULL;
        size_t len = 0;

        label i = 0;
        while ((getline(&line, &len, fp)) != -1) {
            // using printf() in all tests for consistency

            sscanf(line, "%lg", &fld[i]);
            i++;
        }
        fclose(fp);
        if (line)
            free(line);
    }
    Pout<< "Read C " << fld.size() << " in " << timer.cpuTimeIncrement()
        << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
