/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Test-objectRegistry

Description
    Simple test of objectRegistry functionality.
    Particular focus on the behaviour of subRegistry.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IOstreams.H"
#include "FlatOutput.H"
#include "objectRegistry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// file variable, needed for switching the default in lookupObject etc.
bool recursive = false;


void printRegistry
(
    Foam::Ostream& os,
    const Foam::objectRegistry& obr,
    Foam::label indent = 4
);


void printRegistry
(
    Foam::Ostream& os,
    const Foam::objectRegistry& obr,
    Foam::label indent
)
{
    wordList names(obr.sortedNames());
    wordList regs(obr.sortedNames<objectRegistry>());

    std::string prefix;
    for (label i=indent; i; --i)
    {
        prefix += ' ';
    }

    os  << '#' << prefix.c_str() << obr.name()
        << " parent:" << obr.parent().name() << nl;

    os  << ' ' << prefix.c_str() << "objects: " << flatOutput(names) << nl;
    os  << ' ' << prefix.c_str() << "registries: " << flatOutput(regs) << nl;


    // Print, but skip expansion of sub-registries for now
    for (const word& name : names)
    {
        os  << (regs.found(name) ? '-' : ' ')
            << prefix.c_str() << name << " => " << obr[name]->type() << nl;
    }
    for (label i=indent; i; --i)
    {
        os  << '-'; // divider
    }
    os  << '\n';

    // Now descend into the sub-registries
    for (const word& name : regs)
    {
        const objectRegistry& next = obr.lookupObject<objectRegistry>
        (
            name,
            recursive
        );

        os  << prefix.c_str()
            << "current:" << obr.name() << " next:"
            << next.name() << " next-parent:" << next.parent().name() << nl;

        os  << prefix.c_str() << name << " => " << obr[name]->type();

        if ("dictionary" == obr[name]->type())
        {
            os  << " (skip dictionary)" << nl;
        }
        else
        {
            os  << nl;
            printRegistry(os, next, indent + 4);
        }
    }
}

//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addBoolOption
    (
        "mesh",
        "test with polyMesh objectRegistry instead of runTime"
    );
    argList::addBoolOption
    (
        "skip",
        "skip some parts"
    );
    argList::addArgument("recursive (true|false)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    recursive = Switch(args[1]);

    const bool optMesh   = args.found("mesh");
    const bool optSkip   = args.found("skip");
    const objectRegistry& db = (optMesh ? mesh.thisDb() : runTime);

    Info<<"## start ##" << nl;
    Info<< "db isTime:" << db.isTimeDb() << " addr:" << Foam::name(&db)
        << " time: " << Foam::name(&(db.time())) << nl;

    printRegistry(Info, db);
    Info<< nl;

    const label nRegs = 3;

    // Add some items
    for (label j = 0; j < 3; ++j)
    {
        word entryName = "entry" + name(j);
        db.subRegistry
        (
            entryName,
            true,
            recursive
        );
    }

    Info<<"## initially populated ##" << nl;
    printRegistry(Info, db);
    Info<< nl;


    // create a few sub-registries
    for (label i = 0; i < nRegs; ++i)
    {
        word regName = "subreg" + name(i);

        const objectRegistry& subreg = db.subRegistry
        (
            regName,
            true,
            recursive
        );

        for (label j = 0; j < 3; ++j)
        {
            word entryName  = "entry" + name(j);

            subreg.subRegistry
            (
                entryName,
                true,
                recursive
            );
            subreg.subRegistry
            (
                "$" + entryName, // qualified to avoid collisions
                true,
                recursive
            );
        }
    }

    Info<<"## after adding sub-registries" << nl;
    printRegistry(Info, db);
    Info<< nl;

    // Add further items into top-level
    for (label j = 0; j < 6; ++j)
    {
        word entryName = "entry" + name(j);
        db.subRegistry
        (
            entryName,
            true,
            recursive
        );
    }

    Info<< "after adding some entries, top-level now contains: "
        << flatOutput(db.names()) << endl;

    Info<<"## Now attempt to add a few more entries ##" << nl;

    // Try adding the same items into sub registry
    // create a few sub-registries
    for (label i = 0; i < nRegs; ++i)
    {
        word regName = "subreg" + name(i);

        const objectRegistry& subreg = db.subRegistry
        (
            regName,
            false,
            recursive
        );

        if (!optSkip)
        {
            for (label j = 0; j < 6; ++j)
            {
                word entryName = "entry" + name(j);

                subreg.subRegistry
                (
                    entryName,
                    true,
                    recursive
                );
            }
        }
    }

    Info<<"## Complete picture ##" << nl;
    printRegistry(Info, db);
    Info<< nl;

    return 0;
}


// ************************************************************************* //
