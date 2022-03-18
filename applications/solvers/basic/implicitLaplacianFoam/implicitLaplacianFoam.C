/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    laplacianFoam

Group
    grpBasicSolvers

Description
    Explict version of laplacian foam

Description from Laplacian Foam
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "Random.H"
#include "clock.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "postProcess.H"

//    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;


    // use max value to reduce operation cound
    dimensionedScalar DT_max("DT0", DT.dimensions(), scalar(1.0));

   // Random rnd(clock::getTime());
    //Info<<"Random = "<<rnd.scalar01()<<endl;

        Info << "  ClockTime before = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    //while (simple.loop(runTime))  //check simple.loop inside
       
    volScalarField Tinit=T;

    while (runTime.loop())
    {
/*
        while (simple.correctNonOrthogonal())
        {
            
            fvScalarMatrix TEqn
            (
                //fvm::ddt(T) - fvc::laplacian(DT, T)  // reduce operation count
                fvm::ddt(T) - DT_max*fvc::laplacian(T)
                //==
                //fvOptions(T)                         // reduce operation count
            );

            //fvOptions.constrain(TEqn); // reduce operation count
            TEqn.solve();
            //fvOptions.correct(T); // reduce operation count


            // eliminate (maybe) construction of fvScalarMatrix
            solve(fvm::ddt(T) - DT_max*fvc::laplacian(T));
        }
*/
        //#include "write.H"

        //runTime.printExecutionTime(Info);
//        runTime.write();

        solve(fvm::ddt(T) - DT_max*fvm::laplacian(T));
        T=Tinit;
     //   T=Tinit*(1+0.05*rnd.scalar01());

      //  Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"

    }

          Info << "  ClockTime after = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
