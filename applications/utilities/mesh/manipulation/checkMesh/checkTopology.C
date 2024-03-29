/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "checkTopology.H"
#include "polyMesh.H"
#include "Time.H"
#include "regionSplit.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOmanip.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "vtkSurfaceWriter.H"
#include "checkTools.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class PatchType>
void Foam::checkPatch
(
    const bool allGeometry,
    const word& name,
    const PatchType& pp,
    pointSet& points
)
{
    Info<< "    "
        << setw(20) << name
        << setw(9) << returnReduce(pp.size(), sumOp<label>())
        << setw(9) << returnReduce(pp.nPoints(), sumOp<label>());

    if (!Pstream::parRun())
    {
        typedef typename PatchType::surfaceTopo TopoType;
        TopoType pTyp = pp.surfaceType();

        if (pp.empty())
        {
            Info<< setw(34) << "ok (empty)";
        }
        else if (pTyp == TopoType::MANIFOLD)
        {
            if (pp.checkPointManifold(true, &points))
            {
                Info<< setw(34)
                    << "multiply connected (shared point)";
            }
            else
            {
                Info<< setw(34) << "ok (closed singly connected)";
            }

            // Add points on non-manifold edges to make set complete
            pp.checkTopology(false, &points);
        }
        else
        {
            pp.checkTopology(false, &points);

            if (pTyp == TopoType::OPEN)
            {
                Info<< setw(34)
                    << "ok (non-closed singly connected)";
            }
            else
            {
                Info<< setw(34)
                    << "multiply connected (shared edge)";
            }
        }
    }

    if (allGeometry)
    {
        const labelList& mp = pp.meshPoints();

        if (returnReduce(mp.size(), sumOp<label>()) > 0)
        {
            boundBox bb(pp.points(), mp, true); // reduce
            Info<< ' ' << bb;
        }
    }
}


Foam::label Foam::checkTopology
(
    const polyMesh& mesh,
    const bool allTopology,
    const bool allGeometry,
    autoPtr<surfaceWriter>& surfWriter,
    const autoPtr<writer<scalar>>& setWriter
)
{
    label noFailedChecks = 0;

    Info<< "Checking topology..." << endl;

    // Check if the boundary definition is unique
    mesh.boundaryMesh().checkDefinition(true);

    // Check that empty patches cover all sides of the mesh
    {
        label nEmpty = 0;
        forAll(mesh.boundaryMesh(), patchi)
        {
            if (isA<emptyPolyPatch>(mesh.boundaryMesh()[patchi]))
            {
                nEmpty += mesh.boundaryMesh()[patchi].size();
            }
        }
        reduce(nEmpty, sumOp<label>());
        label nTotCells = returnReduce(mesh.cells().size(), sumOp<label>());

        // These are actually warnings, not errors.
        if (nTotCells && (nEmpty % nTotCells))
        {
            Info<< " ***Total number of faces on empty patches"
                << " is not divisible by the number of cells in the mesh."
                << " Hence this mesh is not 1D or 2D."
                << endl;
        }
    }

    // Check if the boundary processor patches are correct
    mesh.boundaryMesh().checkParallelSync(true);

    // Check names of zones are equal
    mesh.cellZones().checkDefinition(true);
    if (mesh.cellZones().checkParallelSync(true))
    {
        noFailedChecks++;
    }
    mesh.faceZones().checkDefinition(true);
    if (mesh.faceZones().checkParallelSync(true))
    {
        noFailedChecks++;
    }
    mesh.pointZones().checkDefinition(true);
    if (mesh.pointZones().checkParallelSync(true))
    {
        noFailedChecks++;
    }


    {
        cellSet cells(mesh, "illegalCells", mesh.nCells()/100);

        forAll(mesh.cells(), celli)
        {
            const cell& cFaces = mesh.cells()[celli];

            if (cFaces.size() <= 3)
            {
                cells.insert(celli);
            }
            for (const label facei : cFaces)
            {
                if (facei < 0 || facei >= mesh.nFaces())
                {
                    cells.insert(celli);
                    break;
                }
            }
        }
        const label nCells = returnReduce(cells.size(), sumOp<label>());

        if (nCells > 0)
        {
            Info<< "    Illegal cells (less than 4 faces or out of range faces)"
                << " found,  number of cells: " << nCells << endl;
            noFailedChecks++;

            Info<< "  <<Writing " << nCells
                << " illegal cells to set " << cells.name() << endl;
            cells.instance() = mesh.pointsInstance();
            cells.write();
            if (surfWriter)
            {
                mergeAndWrite(*surfWriter, cells);
            }
        }
        else
        {
            Info<< "    Cell to face addressing OK." << endl;
        }
    }


    {
        pointSet points(mesh, "unusedPoints", mesh.nPoints()/100);
        if (mesh.checkPoints(true, &points))
        {
            noFailedChecks++;

            label nPoints = returnReduce(points.size(), sumOp<label>());

            Info<< "  <<Writing " << nPoints
                << " unused points to set " << points.name() << endl;
            points.instance() = mesh.pointsInstance();
            points.write();
            if (setWriter)
            {
                mergeAndWrite(*setWriter, points);
            }
        }
    }

    {
        faceSet faces(mesh, "upperTriangularFace", mesh.nFaces()/100);
        if (mesh.checkUpperTriangular(true, &faces))
        {
            noFailedChecks++;
        }

        const label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            Info<< "  <<Writing " << nFaces
                << " unordered faces to set " << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (surfWriter)
            {
                mergeAndWrite(*surfWriter, faces);
            }
        }
    }

    {
        faceSet faces(mesh, "outOfRangeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceVertices(true, &faces))
        {
            noFailedChecks++;

            const label nFaces = returnReduce(faces.size(), sumOp<label>());

            Info<< "  <<Writing " << nFaces
                << " faces with out-of-range or duplicate vertices to set "
                << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (surfWriter)
            {
                mergeAndWrite(*surfWriter, faces);
            }
        }
    }

    if (allTopology)
    {
        cellSet cells(mesh, "zipUpCells", mesh.nCells()/100);
        if (mesh.checkCellsZipUp(true, &cells))
        {
            noFailedChecks++;

            const label nCells = returnReduce(cells.size(), sumOp<label>());

            Info<< "  <<Writing " << nCells
                << " cells with over used edges to set " << cells.name()
                << endl;
            cells.instance() = mesh.pointsInstance();
            cells.write();
            if (surfWriter)
            {
                mergeAndWrite(*surfWriter, cells);
            }

        }
    }

    if (allTopology)
    {
        faceSet faces(mesh, "edgeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceFaces(true, &faces))
        {
            noFailedChecks++;
        }

        const label nFaces = returnReduce(faces.size(), sumOp<label>());
        if (nFaces > 0)
        {
            Info<< "  <<Writing " << nFaces
                << " faces with non-standard edge connectivity to set "
                << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (surfWriter)
            {
                mergeAndWrite(*surfWriter, faces);
            }
        }
    }

    if (allTopology)
    {
        labelList nInternalFaces(mesh.nCells(), Zero);

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            nInternalFaces[mesh.faceOwner()[facei]]++;
            nInternalFaces[mesh.faceNeighbour()[facei]]++;
        }
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        forAll(patches, patchi)
        {
            if (patches[patchi].coupled())
            {
                const labelUList& owners = patches[patchi].faceCells();

                for (const label facei : owners)
                {
                    nInternalFaces[facei]++;
                }
            }
        }

        cellSet oneCells(mesh, "oneInternalFaceCells", mesh.nCells()/100);
        cellSet twoCells(mesh, "twoInternalFacesCells", mesh.nCells()/100);

        forAll(nInternalFaces, celli)
        {
            if (nInternalFaces[celli] <= 1)
            {
                oneCells.insert(celli);
            }
            else if (nInternalFaces[celli] == 2)
            {
                twoCells.insert(celli);
            }
        }

        label nOneCells = returnReduce(oneCells.size(), sumOp<label>());

        if (nOneCells > 0)
        {
            Info<< "  <<Writing " << nOneCells
                << " cells with zero or one non-boundary face to set "
                << oneCells.name()
                << endl;
            oneCells.instance() = mesh.pointsInstance();
            oneCells.write();
            if (surfWriter)
            {
                mergeAndWrite(*surfWriter, oneCells);
            }
        }

        label nTwoCells = returnReduce(twoCells.size(), sumOp<label>());

        if (nTwoCells > 0)
        {
            Info<< "  <<Writing " << nTwoCells
                << " cells with two non-boundary faces to set "
                << twoCells.name()
                << endl;
            twoCells.instance() = mesh.pointsInstance();
            twoCells.write();
            if (surfWriter)
            {
                mergeAndWrite(*surfWriter, twoCells);
            }
        }
    }

    {
        regionSplit rs(mesh);

        if (rs.nRegions() <= 1)
        {
            Info<< "    Number of regions: " << rs.nRegions() << " (OK)."
                << endl;

        }
        else
        {
            Info<< "   *Number of regions: " << rs.nRegions() << endl;

            Info<< "    The mesh has multiple regions which are not connected "
                   "by any face." << endl
                << "  <<Writing region information to "
                << mesh.time().timeName()/"cellToRegion"
                << endl;

            labelIOList ctr
            (
                IOobject
                (
                    "cellToRegion",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rs
            );
            ctr.write();


            // Points in multiple regions
            pointSet points
            (
                mesh,
                "multiRegionPoints",
                mesh.nPoints()/1000
            );

            // Is region disconnected
            boolList regionDisconnected(rs.nRegions(), true);
            if (allTopology)
            {
                // -1   : not assigned
                // -2   : multiple regions
                // >= 0 : single region
                labelList pointToRegion(mesh.nPoints(), -1);

                for
                (
                    label facei = mesh.nInternalFaces();
                    facei < mesh.nFaces();
                    ++facei
                )
                {
                    const label regioni = rs[mesh.faceOwner()[facei]];
                    const face& f = mesh.faces()[facei];
                    for (const label verti : f)
                    {
                        label& pRegion = pointToRegion[verti];
                        if (pRegion == -1)
                        {
                            pRegion = regioni;
                        }
                        else if (pRegion == -2)
                        {
                            // Already marked
                            regionDisconnected[regioni] = false;
                        }
                        else if (pRegion != regioni)
                        {
                            // Multiple regions
                            regionDisconnected[regioni] = false;
                            regionDisconnected[pRegion] = false;
                            pRegion = -2;
                            points.insert(verti);
                        }
                    }
                }

                Pstream::listCombineGather(regionDisconnected, andEqOp<bool>());
                Pstream::listCombineScatter(regionDisconnected);
            }



            // write cellSet for each region
            PtrList<cellSet> cellRegions(rs.nRegions());
            for (label i = 0; i < rs.nRegions(); i++)
            {
                cellRegions.set
                (
                    i,
                    new cellSet
                    (
                        mesh,
                        "region" + Foam::name(i),
                        mesh.nCells()/100
                    )
                );
            }

            forAll(rs, i)
            {
                cellRegions[rs[i]].insert(i);
            }

            for (label i = 0; i < rs.nRegions(); i++)
            {
                Info<< "  <<Writing region " << i;
                if (allTopology)
                {
                    if (regionDisconnected[i])
                    {
                        Info<< " (fully disconnected)";
                    }
                    else
                    {
                        Info<< " (point connected)";
                    }
                }
                Info<< " with "
                    << returnReduce(cellRegions[i].size(), sumOp<scalar>())
                    << " cells to cellSet " << cellRegions[i].name() << endl;

                cellRegions[i].write();
            }

            label nPoints = returnReduce(points.size(), sumOp<label>());
            if (nPoints)
            {
                Info<< "  <<Writing " << nPoints
                    << " points that are in multiple regions to set "
                    << points.name() << endl;
                points.write();
                if (setWriter)
                {
                    mergeAndWrite(*setWriter, points);
                }
            }
        }
    }

    // Non-manifold points
    pointSet points
    (
        mesh,
        "nonManifoldPoints",
        mesh.nPoints()/1000
    );

    {
        if (!Pstream::parRun())
        {
            Info<< "\nChecking patch topology for multiply connected"
                << " surfaces..." << endl;
        }
        else
        {
            Info<< "\nChecking basic patch addressing..." << endl;
        }

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        Pout.setf(ios_base::left);

        Info<< "    "
            << setw(20) << "Patch"
            << setw(9) << "Faces"
            << setw(9) << "Points";
        if (!Pstream::parRun())
        {
            Info<< setw(34) << "Surface topology";
        }
        if (allGeometry)
        {
            Info<< " Bounding box";
        }
        Info<< endl;

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (!isA<processorPolyPatch>(pp))
            {
                checkPatch(allGeometry, pp.name(), pp, points);
                Info<< endl;
            }
        }

        //Info.setf(ios_base::right);
    }

    {
        if (!Pstream::parRun())
        {
            Info<< "\nChecking faceZone topology for multiply connected"
                << " surfaces..." << endl;
        }
        else
        {
            Info<< "\nChecking basic faceZone addressing..." << endl;
        }

        Pout.setf(ios_base::left);

        const faceZoneMesh& faceZones = mesh.faceZones();

        if (faceZones.size())
        {
            Info<< "    "
                << setw(20) << "FaceZone"
                << setw(9) << "Faces"
                << setw(9) << "Points";

            if (!Pstream::parRun())
            {
                Info<< setw(34) << "Surface topology";
            }
            if (allGeometry)
            {
                Info<< " Bounding box";
            }
            Info<< endl;

            for (const faceZone& fz : faceZones)
            {
                checkPatch(allGeometry, fz.name(), fz(), points);
                Info<< endl;
            }
        }
        else
        {
            Info<< "    No faceZones found."<<endl;
        }
    }

    const label nPoints = returnReduce(points.size(), sumOp<label>());

    if (nPoints)
    {
        Info<< "  <<Writing " << nPoints
            << " conflicting points to set " << points.name() << endl;
        points.instance() = mesh.pointsInstance();
        points.write();
        if (setWriter)
        {
            mergeAndWrite(*setWriter, points);
        }
    }

    {
        Info<< "\nChecking basic cellZone addressing..." << endl;

        Pout.setf(ios_base::left);

        const cellZoneMesh& cellZones = mesh.cellZones();

        if (cellZones.size())
        {
            Info<< "    "
                << setw(20) << "CellZone"
                << setw(13) << "Cells"
                << setw(13) << "Points"
                << setw(13) << "Volume"
                << "BoundingBox" << endl;

            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const scalarField& cellVolumes = mesh.cellVolumes();

            bitSet isZonePoint(mesh.nPoints());

            for (const cellZone& cZone : cellZones)
            {
                boundBox bb;
                isZonePoint.reset();  // clears all bits (reset count)
                scalar v = 0.0;

                for (const label celli : cZone)
                {
                    v += cellVolumes[celli];
                    for (const label facei : cells[celli])
                    {
                        const face& f = faces[facei];
                        for (const label verti : f)
                        {
                            if (isZonePoint.set(verti))
                            {
                                bb.add(mesh.points()[verti]);
                            }
                         }
                    }
                }

                bb.reduce();  // Global min/max

                Info<< "    "
                    << setw(19) << cZone.name()
                    << ' ' << setw(12)
                    << returnReduce(cZone.size(), sumOp<label>())
                    << ' ' << setw(12)
                    << returnReduce(isZonePoint.count(), sumOp<label>())
                    << ' ' << setw(12)
                    << returnReduce(v, sumOp<scalar>())
                    << ' ' << bb << endl;
            }
        }
        else
        {
            Info<< "    No cellZones found."<<endl;
        }
    }

    // Force creation of all addressing if requested.
    // Errors will be reported as required
    if (allTopology)
    {
        mesh.cells();
        mesh.faces();
        mesh.edges();
        mesh.points();
        mesh.faceOwner();
        mesh.faceNeighbour();
        mesh.cellCells();
        mesh.edgeCells();
        mesh.pointCells();
        mesh.edgeFaces();
        mesh.pointFaces();
        mesh.cellEdges();
        mesh.faceEdges();
        mesh.pointEdges();
    }

    return noFailedChecks;
}


// ************************************************************************* //
