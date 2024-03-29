/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "UPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::commsStruct::commsStruct()
:
    above_(-1),
    below_(),
    allBelow_(),
    allNotBelow_()
{}


Foam::UPstream::commsStruct::commsStruct
(
    const label above,
    const labelList& below,
    const labelList& allBelow,
    const labelList& allNotBelow
)
:
    above_(above),
    below_(below),
    allBelow_(allBelow),
    allNotBelow_(allNotBelow)
{}


Foam::UPstream::commsStruct::commsStruct
(
    const label nProcs,
    const label myProcID,
    const label above,
    const labelList& below,
    const labelList& allBelow
)
:
    above_(above),
    below_(below),
    allBelow_(allBelow),
    allNotBelow_(nProcs - allBelow.size() - 1)
{
    boolList inBelow(nProcs, false);

    forAll(allBelow, belowI)
    {
        inBelow[allBelow[belowI]] = true;
    }

    label notI = 0;
    forAll(inBelow, proci)
    {
        if ((proci != myProcID) && !inBelow[proci])
        {
            allNotBelow_[notI++] = proci;
        }
    }
    if (notI != allNotBelow_.size())
    {
        FatalErrorInFunction << "problem!" << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::UPstream::commsStruct::operator==(const commsStruct& comm) const
{
    return
    (
        (above_ == comm.above())
     && (below_ == comm.below())
     && (allBelow_ == allBelow())
     && (allNotBelow_ == allNotBelow())
    );
}


bool Foam::UPstream::commsStruct::operator!=(const commsStruct& comm) const
{
    return !operator==(comm);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const UPstream::commsStruct& comm)
{
    os  << comm.above_ << token::SPACE
        << comm.below_ << token::SPACE
        << comm.allBelow_ << token::SPACE
        << comm.allNotBelow_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
