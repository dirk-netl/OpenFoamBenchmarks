/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "boundaryAdjointContribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundaryAdjointContribution, 0);
defineRunTimeSelectionTable(boundaryAdjointContribution, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundaryAdjointContribution::boundaryAdjointContribution
(
    const word& managerName,
    const word& adjointSolverName,
    const word& simulationType,
    const fvPatch& patch
)
:
    patch_(patch)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<boundaryAdjointContribution> boundaryAdjointContribution::New
(
    const word& managerName,
    const word& adjointSolverName,
    const word& simulationType,
    const fvPatch& patch
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(simulationType);

    if (!cstrIter.found())
    {
        FatalErrorInLookup
        (
            "boundaryAdjointContribution",
            simulationType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalError);
    }

    return
        autoPtr<boundaryAdjointContribution>
        (
            cstrIter()
            (
                managerName,
                adjointSolverName,
                simulationType,
                patch
            )
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> boundaryAdjointContribution::adjointTMVariable1Source()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::adjointTMVariable2Source()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::TMVariable1Diffusion()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::TMVariable2Diffusion()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::TMVariable1()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


tmp<scalarField> boundaryAdjointContribution::TMVariable2()
{
    return tmp<scalarField>::New(patch_.size(), Zero);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //