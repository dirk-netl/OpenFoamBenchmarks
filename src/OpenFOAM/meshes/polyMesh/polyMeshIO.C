/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "Time.H"
#include "cellIOList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::setInstance
(
    const fileName& inst,
    const IOobject::writeOption wOpt
)
{
    DebugInFunction << "Resetting file instance to " << inst << endl;

    points_.writeOpt(wOpt);
    points_.instance() = inst;

    faces_.writeOpt(wOpt);
    faces_.instance() = inst;

    owner_.writeOpt(wOpt);
    owner_.instance() = inst;

    neighbour_.writeOpt(wOpt);
    neighbour_.instance() = inst;

    boundary_.writeOpt(wOpt);
    boundary_.instance() = inst;

    pointZones_.writeOpt(wOpt);
    pointZones_.instance() = inst;

    faceZones_.writeOpt(wOpt);
    faceZones_.instance() = inst;

    cellZones_.writeOpt(wOpt);
    cellZones_.instance() = inst;

    if (tetBasePtIsPtr_)
    {
        tetBasePtIsPtr_->writeOpt(wOpt);
        tetBasePtIsPtr_->instance() = inst;
    }
}


Foam::polyMesh::readUpdateState Foam::polyMesh::readUpdate()
{
    DebugInFunction << "Updating mesh based on saved data." << endl;

    // Find the point and cell instance
    fileName pointsInst(time().findInstance(meshDir(), "points"));
    fileName facesInst(time().findInstance(meshDir(), "faces"));
    //fileName boundaryInst(time().findInstance(meshDir(), "boundary"));

    if (debug)
    {
        Info<< "Faces instance: old = " << facesInstance()
            << " new = " << facesInst << nl
            << "Points instance: old = " << pointsInstance()
            << " new = " << pointsInst << endl;
    }

    if (facesInst != facesInstance())
    {
        // Topological change
        if (debug)
        {
            Info<< "Topological change" << endl;
        }

        clearOut();

        // Set instance to new instance. Note that points instance can differ
        // from from faces instance.
        setInstance(facesInst);
        points_.instance() = pointsInst;

        points_ = pointIOField
        (
            IOobject
            (
                "points",
                pointsInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        faces_ = faceCompactIOList
        (
            IOobject
            (
                "faces",
                facesInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        owner_ = labelIOList
        (
            IOobject
            (
                "owner",
                facesInst,
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        neighbour_ = labelIOList
        (
            IOobject
            (
                "neighbour",
                facesInst,
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        // Reset the boundary patches
        polyBoundaryMesh newBoundary
        (
            IOobject
            (
                "boundary",
                facesInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        // Check that patch types and names are unchanged
        bool boundaryChanged = false;

        if (newBoundary.size() != boundary_.size())
        {
            boundaryChanged = true;
        }
        else
        {
            wordList newTypes = newBoundary.types();
            wordList newNames = newBoundary.names();

            wordList oldTypes = boundary_.types();
            wordList oldNames = boundary_.names();

            forAll(oldTypes, patchi)
            {
                if
                (
                    oldTypes[patchi] != newTypes[patchi]
                 || oldNames[patchi] != newNames[patchi]
                )
                {
                    boundaryChanged = true;
                    break;
                }
            }
        }

        if (boundaryChanged)
        {
            WarningInFunction
                << "Number of patches has changed.  This may have "
                << "unexpected consequences.  Proceed with care." << endl;

            boundary_.clear();
            boundary_.setSize(newBoundary.size());

            forAll(newBoundary, patchi)
            {
                boundary_.set(patchi, newBoundary[patchi].clone(boundary_));
            }
        }
        else
        {
            forAll(boundary_, patchi)
            {
                boundary_[patchi] = polyPatch
                (
                    newBoundary[patchi].name(),
                    newBoundary[patchi].size(),
                    newBoundary[patchi].start(),
                    patchi,
                    boundary_,
                    newBoundary[patchi].physicalType(),
                    newBoundary[patchi].inGroups()
                );
            }
        }


        // Boundary is set so can use initMesh now (uses boundary_ to
        // determine internal and active faces)

        if (!owner_.headerClassName().empty())
        {
            initMesh();
        }
        else
        {
            cellCompactIOList cells
            (
                IOobject
                (
                    "cells",
                    facesInst,
                    meshSubDir,
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

            // Recalculate the owner/neighbour addressing and reset the
            // primitiveMesh
            initMesh(cells);
        }


        // Even if number of patches stayed same still recalculate boundary
        // data.

        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();

        // Derived info
        bounds_ = boundBox(points_);
        geometricD_ = Zero;
        solutionD_ = Zero;

        // Zones
        pointZoneMesh newPointZones
        (
            IOobject
            (
                "pointZones",
                facesInst,
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        label oldSize = pointZones_.size();

        if (newPointZones.size() <= pointZones_.size())
        {
            pointZones_.setSize(newPointZones.size());
        }

        // Reset existing ones
        forAll(pointZones_, czI)
        {
            pointZones_[czI] = newPointZones[czI];
        }

        // Extend with extra ones
        pointZones_.setSize(newPointZones.size());

        for (label czI = oldSize; czI < newPointZones.size(); czI++)
        {
            pointZones_.set(czI, newPointZones[czI].clone(pointZones_));
        }


        faceZoneMesh newFaceZones
        (
            IOobject
            (
                "faceZones",
                facesInst,
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        oldSize = faceZones_.size();

        if (newFaceZones.size() <= faceZones_.size())
        {
            faceZones_.setSize(newFaceZones.size());
        }

        // Reset existing ones
        forAll(faceZones_, fzI)
        {
            faceZones_[fzI].resetAddressing
            (
                newFaceZones[fzI],
                newFaceZones[fzI].flipMap()
            );
        }

        // Extend with extra ones
        faceZones_.setSize(newFaceZones.size());

        for (label fzI = oldSize; fzI < newFaceZones.size(); fzI++)
        {
            faceZones_.set(fzI, newFaceZones[fzI].clone(faceZones_));
        }


        cellZoneMesh newCellZones
        (
            IOobject
            (
                "cellZones",
                facesInst,
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        oldSize = cellZones_.size();

        if (newCellZones.size() <= cellZones_.size())
        {
            cellZones_.setSize(newCellZones.size());
        }

        // Reset existing ones
        forAll(cellZones_, czI)
        {
            cellZones_[czI] = newCellZones[czI];
        }

        // Extend with extra ones
        cellZones_.setSize(newCellZones.size());

        for (label czI = oldSize; czI < newCellZones.size(); czI++)
        {
            cellZones_.set(czI, newCellZones[czI].clone(cellZones_));
        }

        // Re-read tet base points
        tetBasePtIsPtr_ = readTetBasePtIs();


        if (boundaryChanged)
        {
            return polyMesh::TOPO_PATCH_CHANGE;
        }
        else
        {
            return polyMesh::TOPO_CHANGE;
        }
    }
    else if (pointsInst != pointsInstance())
    {
        // Points moved
        if (debug)
        {
            Info<< "Point motion" << endl;
        }

        pointIOField newPoints
        (
            IOobject
            (
                "points",
                pointsInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // Re-read tet base points
        autoPtr<labelIOList> newTetBasePtIsPtr = readTetBasePtIs();

        // Update all geometry
        updateGeomPoints(std::move(newPoints), newTetBasePtIsPtr);

        return polyMesh::POINTS_MOVED;
    }
    else
    {
        if (debug)
        {
            Info<< "No change" << endl;
        }

        return polyMesh::UNCHANGED;
    }
}


// ************************************************************************* //
