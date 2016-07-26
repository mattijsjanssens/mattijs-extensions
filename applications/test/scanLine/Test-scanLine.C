/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 M. Janssens
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
    Test-scanLine

Description
    Scanline algorithm

\*---------------------------------------------------------------------------*/

#include "argList.H"
//#include "Time.H"
#include "OBJstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Line
(
    scalar x1,
    scalar y1,
    scalar x2,
    scalar y2
)
{
    // https://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm#C.2B.2B

    OBJstream str("line.obj");

    // Bresenham's line algorithm
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if(steep)
    {
        Swap(x1, y1);
        Swap(x2, y2);
    }

    if(x1 > x2)
    {
        Swap(x1, x2);
        Swap(y1, y2);
    }
 
    const scalar dx = x2 - x1;
    const scalar dy = fabs(y2 - y1);
 
    scalar error = dx / 2.0f;
    const label ystep = (y1 < y2) ? 1 : -1;
    label y = label(y1);
 
    const label maxX = label(x2);
 
    for(label x=label(x1); x<maxX; x++)
    {
        if (steep)
        {
            //SetPixel(y,x, color);
            str.write(point(y, x, 0));
        }
    	else
        {
            //SetPixel(x,y, color);
            str.write(point(x, y, 0));
        }
 
        error -= dy;
        if (error < 0)
        {
            y += ystep;
            error += dx;
        }
    }
}


int main(int argc, char *argv[])
{
//    #include "setRootCase.H"
//    #include "createTime.H"

    
    Line(0, 0, 10, 10);

    Info<< "end" << endl;
}


// ************************************************************************* //
