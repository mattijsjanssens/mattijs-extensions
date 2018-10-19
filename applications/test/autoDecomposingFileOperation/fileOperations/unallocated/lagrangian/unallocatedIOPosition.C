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

\*---------------------------------------------------------------------------*/

#include "unallocatedIOPosition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unallocatedIOPosition::unallocatedIOPosition(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<unallocatedIOPosition>();

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // Make sure to not check header class name
        readStream(word::null) >> *this;
        close();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::unallocatedIOPosition::write(const bool valid) const
{
    return regIOobject::write(IDLList<basicParticle>::size());
}


bool Foam::unallocatedIOPosition::writeData(Ostream& os) const
{
    os  << IDLList<basicParticle>::size() << nl << token::BEGIN_LIST << nl;

    forAllConstIter(typename IDLList<basicParticle>, *this, iter)
    {
        iter().writePosition(os);
        os  << nl;
    }

    os  << token::END_LIST << endl;

    return os.good();
}


bool Foam::unallocatedIOPosition::readData(Istream& is)
{
    token firstToken(is);

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        is.readBeginList
        (
            "unallocatedIOPosition::readData(Istream&)"
        );

        for (label i=0; i<s; i++)
        {
            // Read position only
            append(new basicParticle(is, false));
        }

        // Read end of contents
        is.readEndList
        (
            "unallocatedIOPosition::readData(Istream&)"
        );
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction
            (
                is
            )   << "incorrect first token, '(', found "
                << firstToken.info() << exit(FatalIOError);
        }

        token lastToken(is);
        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);

            // Read position only
            append(new basicParticle(is, false));
            is  >> lastToken;
        }
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "incorrect first token, expected <int> or '(', found "
            << firstToken.info() << exit(FatalIOError);
    }

    // Check state of IOstream
    is.check
    (
        "void unallocatedIOPosition::readData(Istream&)"
    );

    return is.good();
}


// ************************************************************************* //
