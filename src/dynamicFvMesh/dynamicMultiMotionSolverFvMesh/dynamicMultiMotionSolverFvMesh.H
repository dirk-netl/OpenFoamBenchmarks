/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

Class
    Foam::dynamicMultiMotionSolverFvMesh

Description
    Mesh motion described per cellZone. Individual motion solvers solve
    over whole domain but are only applied per cellZone.

SourceFiles
    dynamicMultiMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMultiMotionSolverFvMesh_H
#define dynamicMultiMotionSolverFvMesh_H

#include "dynamicFvMesh.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class dynamicMultiMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMultiMotionSolverFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- The motion control function
        PtrList<motionSolver> motionPtr_;

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;


    // Private Member Functions

        //- No copy construct
        dynamicMultiMotionSolverFvMesh
        (
            const dynamicMultiMotionSolverFvMesh&
        ) = delete;

        //- No copy assignment
        void operator=(const dynamicMultiMotionSolverFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicMultiMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicMultiMotionSolverFvMesh
        (
            const IOobject& io,
            const bool doInit=true
        );


    //- Destructor
    ~dynamicMultiMotionSolverFvMesh() = default;


    // Member Functions

        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
