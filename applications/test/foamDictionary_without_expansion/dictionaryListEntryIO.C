/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "dictionaryListEntry.H"
#include "keyType.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::dictionaryListEntry::realSize(const dictionary& dict)
{
    if (dict.size() < 1 || dict.first()->keyword() != "FoamFile")
    {
        return dict.size();
    }
    else
    {
        return dict.size() - 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryListEntry::dictionaryListEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    dictionaryEntry
    (
        Foam::name(realSize(parentDict)),
        parentDict,
        dictionary::null
    )
{
    token keywordToken(is);
    label s = keywordToken.labelToken();

    is.readBeginList("List");

    for (label i=0; i<s; i++)
    {
        entry::New(*this, is);
    }
    is.readEndList("List");

    is.fatalCheck
    (
        "dictionaryListEntry::dictionaryListEntry"
        "(const dictionary& parentDict, Istream&)"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dictionaryListEntry::write(Ostream& os) const
{
    os  << nl << indent << size()
        << token::SPACE << "// entry " << keyword() << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;

    // Write contents
    dictionary::write(os, false);

    // Write end delimiter
    os << decrIndent << indent << token::END_LIST << nl;

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const dictionaryListEntry&)");
}


// * * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const dictionaryListEntry& de)
{
    de.write(os);
    return os;
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<dictionaryListEntry>& ip
)
{
    const dictionaryListEntry& e = ip.t_;

    os  << "    dictionaryListEntry '" << e.keyword() << "'" << endl;

    return os;
}


// ************************************************************************* //
