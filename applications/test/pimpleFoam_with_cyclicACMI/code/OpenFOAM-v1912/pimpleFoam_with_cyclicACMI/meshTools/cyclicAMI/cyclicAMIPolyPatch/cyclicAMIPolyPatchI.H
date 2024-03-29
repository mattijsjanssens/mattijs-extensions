/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::wordList& Foam::cyclicAMIPolyPatch::neighbPatchNames() const
{
    if (nbrPatchNames_.empty())
    {
        // Try and use patchGroup to find samplePatch and sampleRegion
        const labelList nbrPatchIDs(coupleGroup_.findOtherPatchIDs(*this));

        nbrPatchNames_.setSize(nbrPatchIDs.size());
        forAll(nbrPatchIDs, i)
        {
            nbrPatchNames_[i] = boundaryMesh()[nbrPatchIDs[i]].name();
        }

        // Note: could store nbrPatchIDs_ directly here and bypass
        //       neighbPatchIDs but that does additional checking

        if (nbrPatchNames_.empty())
        {
            FatalErrorInFunction<< "Did not find other patches "
                << " in patchGroup " << coupleGroup_.name()
                << " when searching for neigubour patches for patch " << name()
                << exit(FatalError);
        }
    }
    return nbrPatchNames_;
}


inline const Foam::word& Foam::cyclicAMIPolyPatch::neighbPatchName() const
{
    return neighbPatchNames()[0];
}


inline const Foam::vector& Foam::cyclicAMIPolyPatch::rotationAxis() const
{
    return rotationAxis_;
}


inline const Foam::point& Foam::cyclicAMIPolyPatch::rotationCentre() const
{
    return rotationCentre_;
}


inline const Foam::vector& Foam::cyclicAMIPolyPatch::separationVector() const
{
    return separationVector_;
}


// ************************************************************************* //
