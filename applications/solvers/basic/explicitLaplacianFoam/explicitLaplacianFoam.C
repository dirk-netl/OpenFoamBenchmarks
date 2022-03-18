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
#include <chrono>

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
    dimensionedScalar center("center", dimensionSet(0,0,0,0,0,0,0), 0.4);
    dimensionedScalar coeff("coeff", dimensionSet(0,0,0,0,0,0,0), 0.1);

        //  Info << "  ClockTime before = " << runTime.elapsedClockTime() << " s"
        //    << nl << endl;


    //while (simple.loop(runTime))  //check simple.loop inside

    volScalarField Tneighbor_sum=T*0.0;
    volScalarField T_temp=T*0.0;
    while (runTime.loop())  //check simple.loop inside
    {
// Start measuring time
     auto begin = std::chrono::high_resolution_clock::now();

 	const fvPatchList& patches = mesh.boundary();
	Tneighbor_sum = 0.0;  
	forAll(T, celli)
	{
 	//Info << patchI<<", "<<T.boundaryField()[patchI].patchInternalField()<<endl;
       
   	const labelList& neighbor = mesh.cellCells()[celli];
        

	forAll(neighbor,i)
 	{
	
	Tneighbor_sum[celli] += T[neighbor[i]]; //Tneighbor_sum += 64bit + cast(T[neighbor[i])
 	}
	}	
        T_temp = center*T;

        T.correctBoundaryConditions();

	forAll(patches, patchI)
	{
	//Info<<patchI<<", "<< T.boundaryField()[patchI].patchInternalField()<<endl;
	forAll(T.boundaryField()[patchI], iFaceLocal)
	{
    	label iFaceGlobal = iFaceLocal + mesh.boundaryMesh()[patchI].start();
    	label adjacentCell = mesh.owner()[iFaceGlobal];
    	Tneighbor_sum[adjacentCell] += T.boundaryField()[patchI][iFaceLocal];
	}
	}
        //solve(fvm::ddt(T) - DT_max*fvc::laplacian(T));
        T = T_temp  + coeff*Tneighbor_sum ; // T=cast(Tneighbor_sum+64bitcenter * cast(T in 64bit))
        //Info<<T<<endl;
       // Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"<<endl;
        
        //comment out to avoid writing data
        //#include "write.H" 
        //runTime.write();
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

          Info << "  Elapsed Time = " << elapsed.count()*1e-9 << " s"
            << nl << endl;
    }


    Info<< "Time = " << runTime.timeName() << nl << endl;
       //   Info << "  ClockTime after = " << runTime.elapsedClockTime() << " s"
       //     << nl << endl;
    //printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
