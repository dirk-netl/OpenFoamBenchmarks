/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a polyMesh for a specified region, or the defaultRegion

Required Variables
    - runTime [Time]

Provided Variables
    - regionName [word]
    - mesh [polyMesh]

\*---------------------------------------------------------------------------*/

Foam::word regionName(Foam::polyMesh::defaultRegion);

{
    if (args.readIfPresent("region", regionName))
    {
        Foam::Info << "Create polyMesh " << regionName;
    }
    else
    {
        Foam::Info << "Create polyMesh";
    }
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;
}


Foam::polyMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

Foam::Info << Foam::endl;

// ************************************************************************* //
