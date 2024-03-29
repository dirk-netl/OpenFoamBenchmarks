/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "domainDecomposition.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "fvMesh.H"
#include "OSspecific.H"
#include "Map.H"
#include "DynamicList.H"
#include "fvFieldDecomposer.H"
#include "IOobjectList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "decompositionModel.H"
#include "hexRef8Data.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::domainDecomposition::mark
(
    const labelList& zoneElems,
    const label zoneI,
    labelList& elementToZone
)
{
    for (const label pointi : zoneElems)
    {
        if (elementToZone[pointi] == -1)
        {
            // First occurrence
            elementToZone[pointi] = zoneI;
        }
        else if (elementToZone[pointi] >= 0)
        {
            // Multiple zones
            elementToZone[pointi] = -2;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::domainDecomposition::domainDecomposition
(
    const IOobject& io,
    const fileName& decompDictFile
)
:
    fvMesh(io),
    facesInstancePointsPtr_
    (
        pointsInstance() != facesInstance()
      ? new pointIOField
        (
            IOobject
            (
                "points",
                facesInstance(),
                polyMesh::meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
      : nullptr
    ),
    decompDictFile_(decompDictFile),
    nProcs_
    (
        decompositionMethod::nDomains
        (
            decompositionModel::New
            (
                *this,
                decompDictFile
            )
        )
    ),
    distributed_(false),
    cellToProc_(nCells()),
    procPointAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    procCellAddressing_(nProcs_),
    procPatchSize_(nProcs_),
    procPatchStartIndex_(nProcs_),
    procNeighbourProcessors_(nProcs_),
    procProcessorPatchSize_(nProcs_),
    procProcessorPatchStartIndex_(nProcs_),
    procProcessorPatchSubPatchIDs_(nProcs_),
    procProcessorPatchSubPatchStarts_(nProcs_)
{
    updateParameters(this->model());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::decompositionModel& Foam::domainDecomposition::model() const
{
    return decompositionModel::New(*this, decompDictFile_);
}


void Foam::domainDecomposition::updateParameters
(
    const dictionary& params
)
{
    params.readIfPresent("distributed", distributed_);
}


bool Foam::domainDecomposition::writeDecomposition(const bool decomposeSets)
{
    Info<< "\nConstructing processor meshes" << endl;

    // Mark point/faces/cells that are in zones.
    // -1   : not in zone
    // -2   : in multiple zones
    // >= 0 : in single given zone
    // This will give direct lookup of elements that are in a single zone
    // and we'll only have to revert back to searching through all zones
    // for the duplicate elements

    // Point zones
    labelList pointToZone(points().size(), -1);

    forAll(pointZones(), zonei)
    {
        mark(pointZones()[zonei], zonei, pointToZone);
    }

    // Face zones
    labelList faceToZone(faces().size(), -1);

    forAll(faceZones(), zonei)
    {
        mark(faceZones()[zonei], zonei, faceToZone);
    }

    // Cell zones
    labelList cellToZone(nCells(), -1);

    forAll(cellZones(), zonei)
    {
        mark(cellZones()[zonei], zonei, cellToZone);
    }


    PtrList<const cellSet> cellSets;
    PtrList<const faceSet> faceSets;
    PtrList<const pointSet> pointSets;
    if (decomposeSets)
    {
        // Read sets
        IOobjectList objects(*this, facesInstance(), "polyMesh/sets");
        {
            IOobjectList cSets(objects.lookupClass(cellSet::typeName));
            forAllConstIters(cSets, iter)
            {
                cellSets.append(new cellSet(*iter()));
            }
        }
        {
            IOobjectList fSets(objects.lookupClass(faceSet::typeName));
            forAllConstIters(fSets, iter)
            {
                faceSets.append(new faceSet(*iter()));
            }
        }
        {
            IOobjectList pSets(objects.lookupClass(pointSet::typeName));
            forAllConstIters(pSets, iter)
            {
                pointSets.append(new pointSet(*iter()));
            }
        }
    }


    // Load refinement data (if any)
    hexRef8Data baseMeshData
    (
        IOobject
        (
            "dummy",
            facesInstance(),
            polyMesh::meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );



    label maxProcCells = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;
    label maxProcFaces = 0;


    // Write out the meshes
    for (label proci = 0; proci < nProcs_; proci++)
    {
        // Create processor points
        const labelList& curPointLabels = procPointAddressing_[proci];

        const pointField& meshPoints = points();

        labelList pointLookup(nPoints(), -1);

        pointField procPoints(curPointLabels.size());

        forAll(curPointLabels, pointi)
        {
            procPoints[pointi] = meshPoints[curPointLabels[pointi]];

            pointLookup[curPointLabels[pointi]] = pointi;
        }

        // Create processor faces
        const labelList& curFaceLabels = procFaceAddressing_[proci];

        const faceList& meshFaces = faces();

        labelList faceLookup(nFaces(), -1);

        faceList procFaces(curFaceLabels.size());

        forAll(curFaceLabels, facei)
        {
            // Mark the original face as used
            // Remember to decrement the index by one (turning index)
            //
            label curF = mag(curFaceLabels[facei]) - 1;

            faceLookup[curF] = facei;

            // get the original face
            labelList origFaceLabels;

            if (curFaceLabels[facei] >= 0)
            {
                // face not turned
                origFaceLabels = meshFaces[curF];
            }
            else
            {
                origFaceLabels = meshFaces[curF].reverseFace();
            }

            // translate face labels into local point list
            face& procFaceLabels = procFaces[facei];

            procFaceLabels.setSize(origFaceLabels.size());

            forAll(origFaceLabels, pointi)
            {
                procFaceLabels[pointi] = pointLookup[origFaceLabels[pointi]];
            }
        }

        // Create processor cells
        const labelList& curCellLabels = procCellAddressing_[proci];

        const cellList& meshCells = cells();

        cellList procCells(curCellLabels.size());

        forAll(curCellLabels, celli)
        {
            const labelList& origCellLabels = meshCells[curCellLabels[celli]];

            cell& curCell = procCells[celli];

            curCell.setSize(origCellLabels.size());

            forAll(origCellLabels, cellFacei)
            {
                curCell[cellFacei] = faceLookup[origCellLabels[cellFacei]];
            }
        }

        // Create processor mesh without a boundary

        fileName processorCasePath
        (
            time().caseName()/("processor" + Foam::name(proci))
        );

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            processorCasePath,
            word("system"),
            word("constant")
        );
        processorDb.setTime(time());

        // create the mesh. Two situations:
        // - points and faces come from the same time ('instance'). The mesh
        //   will get constructed in the same instance.
        // - points come from a different time (moving mesh cases).
        //   It will read the points belonging to the faces instance and
        //   construct the procMesh with it which then gets handled as above.
        //   (so with 'old' geometry).
        //   Only at writing time will it additionally write the current
        //   points.

        autoPtr<polyMesh> procMeshPtr;

        if (facesInstancePointsPtr_)
        {
            // Construct mesh from facesInstance.
            pointField facesInstancePoints
            (
                facesInstancePointsPtr_(),
                curPointLabels
            );

            procMeshPtr = autoPtr<polyMesh>::New
            (
                IOobject
                (
                    this->polyMesh::name(), // region of undecomposed mesh
                    facesInstance(),
                    processorDb,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                std::move(facesInstancePoints),
                std::move(procFaces),
                std::move(procCells)
            );
        }
        else
        {
            procMeshPtr = autoPtr<polyMesh>::New
            (
                IOobject
                (
                    this->polyMesh::name(), // region of undecomposed mesh
                    facesInstance(),
                    processorDb,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                std::move(procPoints),
                std::move(procFaces),
                std::move(procCells)
            );
        }
        polyMesh& procMesh = procMeshPtr();


        // Create processor boundary patches
        const labelList& curPatchSizes = procPatchSize_[proci];

        const labelList& curPatchStarts = procPatchStartIndex_[proci];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[proci];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[proci];

        const labelList& curProcessorPatchStarts =
            procProcessorPatchStartIndex_[proci];

        const labelListList& curSubPatchIDs =
            procProcessorPatchSubPatchIDs_[proci];

        const labelListList& curSubStarts =
            procProcessorPatchSubPatchStarts_[proci];

        const polyPatchList& meshPatches = boundaryMesh();


        // Count the number of inter-proc patches
        label nInterProcPatches = 0;
        forAll(curSubPatchIDs, procPatchi)
        {
            nInterProcPatches += curSubPatchIDs[procPatchi].size();
        }

        PtrList<polyPatch> procPatches
        (
            curPatchSizes.size() + nInterProcPatches
        );

        label nPatches = 0;

        forAll(curPatchSizes, patchi)
        {
            // Get the face labels consistent with the field mapping
            // (reuse the patch field mappers)
            const polyPatch& meshPatch = meshPatches[patchi];

            fvFieldDecomposer::patchFieldDecomposer patchMapper
            (
                SubList<label>
                (
                    curFaceLabels,
                    curPatchSizes[patchi],
                    curPatchStarts[patchi]
                ),
                meshPatch.start()
            );

            // Map existing patches
            procPatches.set
            (
                nPatches,
                meshPatch.clone
                (
                    procMesh.boundaryMesh(),
                    nPatches,
                    patchMapper.directAddressing(),
                    curPatchStarts[patchi]
                )
            );

            nPatches++;
        }

        forAll(curProcessorPatchSizes, procPatchi)
        {
            const labelList& subPatchID = curSubPatchIDs[procPatchi];
            const labelList& subStarts = curSubStarts[procPatchi];

            label curStart = curProcessorPatchStarts[procPatchi];

            forAll(subPatchID, i)
            {
                label size =
                (
                    i < subPatchID.size()-1
                  ? subStarts[i+1] - subStarts[i]
                  : curProcessorPatchSizes[procPatchi] - subStarts[i]
                );

                if (subPatchID[i] == -1)
                {
                    // From internal faces
                    procPatches.set
                    (
                        nPatches,
                        new processorPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi]
                        )
                    );
                }
                else
                {
                    const coupledPolyPatch& pcPatch
                        = refCast<const coupledPolyPatch>
                          (
                              boundaryMesh()[subPatchID[i]]
                          );

                    procPatches.set
                    (
                        nPatches,
                        new processorCyclicPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi],
                            pcPatch.name(),
                            pcPatch.transform()
                        )
                    );
                }

                curStart += size;
                ++nPatches;
            }
        }

        // Add boundary patches
        procMesh.addPatches(procPatches);

        // Create and add zones

        // Point zones
        {
            const pointZoneMesh& pz = pointZones();

            // Go through all the zoned points and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zonePoints(pz.size());

            // Estimate size
            forAll(zonePoints, zonei)
            {
                zonePoints[zonei].setCapacity(pz[zonei].size()/nProcs_);
            }

            // Use the pointToZone map to find out the single zone (if any),
            // use slow search only for shared points.
            forAll(curPointLabels, pointi)
            {
                label curPoint = curPointLabels[pointi];

                label zonei = pointToZone[curPoint];

                if (zonei >= 0)
                {
                    // Single zone.
                    zonePoints[zonei].append(pointi);
                }
                else if (zonei == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(pz, zonei)
                    {
                        label index = pz[zonei].whichPoint(curPoint);

                        if (index != -1)
                        {
                            zonePoints[zonei].append(pointi);
                        }
                    }
                }
            }

            procMesh.pointZones().clearAddressing();
            procMesh.pointZones().setSize(zonePoints.size());
            forAll(zonePoints, zonei)
            {
                procMesh.pointZones().set
                (
                    zonei,
                    pz[zonei].clone
                    (
                        procMesh.pointZones(),
                        zonei,
                        zonePoints[zonei].shrink()
                    )
                );
            }

            if (pz.size())
            {
                // Force writing on all processors
                procMesh.pointZones().writeOpt(IOobject::AUTO_WRITE);
            }
        }

        // Face zones
        {
            const faceZoneMesh& fz = faceZones();

            // Go through all the zoned face and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zoneFaces(fz.size());
            List<DynamicList<bool>> zoneFaceFlips(fz.size());

            // Estimate size
            forAll(zoneFaces, zonei)
            {
                label procSize = fz[zonei].size() / nProcs_;

                zoneFaces[zonei].setCapacity(procSize);
                zoneFaceFlips[zonei].setCapacity(procSize);
            }

            // Go through all the zoned faces and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            forAll(curFaceLabels, facei)
            {
                // Remember to decrement the index by one (turning index)
                //
                label curF = mag(curFaceLabels[facei]) - 1;

                label zonei = faceToZone[curF];

                if (zonei >= 0)
                {
                    // Single zone. Add the face
                    zoneFaces[zonei].append(facei);

                    label index = fz[zonei].whichFace(curF);

                    bool flip = fz[zonei].flipMap()[index];

                    if (curFaceLabels[facei] < 0)
                    {
                        flip = !flip;
                    }

                    zoneFaceFlips[zonei].append(flip);
                }
                else if (zonei == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(fz, zonei)
                    {
                        label index = fz[zonei].whichFace(curF);

                        if (index != -1)
                        {
                            zoneFaces[zonei].append(facei);

                            bool flip = fz[zonei].flipMap()[index];

                            if (curFaceLabels[facei] < 0)
                            {
                                flip = !flip;
                            }

                            zoneFaceFlips[zonei].append(flip);
                        }
                    }
                }
            }

            procMesh.faceZones().clearAddressing();
            procMesh.faceZones().setSize(zoneFaces.size());
            forAll(zoneFaces, zonei)
            {
                procMesh.faceZones().set
                (
                    zonei,
                    fz[zonei].clone
                    (
                        zoneFaces[zonei].shrink(),          // addressing
                        zoneFaceFlips[zonei].shrink(),      // flipmap
                        zonei,
                        procMesh.faceZones()
                    )
                );
            }

            if (fz.size())
            {
                // Force writing on all processors
                procMesh.faceZones().writeOpt(IOobject::AUTO_WRITE);
            }
        }

        // Cell zones
        {
            const cellZoneMesh& cz = cellZones();

            // Go through all the zoned cells and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zoneCells(cz.size());

            // Estimate size
            forAll(zoneCells, zonei)
            {
                zoneCells[zonei].setCapacity(cz[zonei].size()/nProcs_);
            }

            forAll(curCellLabels, celli)
            {
                label curCelli = curCellLabels[celli];

                label zonei = cellToZone[curCelli];

                if (zonei >= 0)
                {
                    // Single zone.
                    zoneCells[zonei].append(celli);
                }
                else if (zonei == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(cz, zonei)
                    {
                        label index = cz[zonei].whichCell(curCelli);

                        if (index != -1)
                        {
                            zoneCells[zonei].append(celli);
                        }
                    }
                }
            }

            procMesh.cellZones().clearAddressing();
            procMesh.cellZones().setSize(zoneCells.size());
            forAll(zoneCells, zonei)
            {
                procMesh.cellZones().set
                (
                    zonei,
                    cz[zonei].clone
                    (
                        zoneCells[zonei].shrink(),
                        zonei,
                        procMesh.cellZones()
                    )
                );
            }

            if (cz.size())
            {
                // Force writing on all processors
                procMesh.cellZones().writeOpt(IOobject::AUTO_WRITE);
            }
        }

        // Set the precision of the points data to be min 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        procMesh.write();

        // Write points if pointsInstance differing from facesInstance
        if (facesInstancePointsPtr_)
        {
            pointIOField pointsInstancePoints
            (
                IOobject
                (
                    "points",
                    pointsInstance(),
                    polyMesh::meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                std::move(procPoints)
            );
            pointsInstancePoints.write();
        }


        // Decompose any sets
        if (decomposeSets)
        {
            forAll(cellSets, i)
            {
                const cellSet& cs = cellSets[i];
                cellSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(curCellLabels, i)
                {
                    if (cs.found(curCellLabels[i]))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
            forAll(faceSets, i)
            {
                const faceSet& cs = faceSets[i];
                faceSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(curFaceLabels, i)
                {
                    if (cs.found(mag(curFaceLabels[i])-1))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
            forAll(pointSets, i)
            {
                const pointSet& cs = pointSets[i];
                pointSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(curPointLabels, i)
                {
                    if (cs.found(curPointLabels[i]))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
        }


        // Optional hexRef8 data
        hexRef8Data
        (
            IOobject
            (
                "dummy",
                facesInstance(),
                polyMesh::meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            baseMeshData,
            procCellAddressing_[proci],
            procPointAddressing_[proci]
        ).write();


        // Statistics
        Info<< nl << "Processor " << proci;

        if (procMesh.nCells())
        {
            Info<< nl << "    ";
        }
        else
        {
            Info<< ": ";
        }

        Info<< "Number of cells = " << procMesh.nCells() << nl;

        maxProcCells = max(maxProcCells, procMesh.nCells());

        label nBoundaryFaces = 0;
        label nProcPatches = 0;
        label nProcFaces = 0;

        forAll(procMesh.boundaryMesh(), patchi)
        {
            if (isA<processorPolyPatch>(procMesh.boundaryMesh()[patchi]))
            {
                const processorPolyPatch& ppp =
                refCast<const processorPolyPatch>
                (
                    procMesh.boundaryMesh()[patchi]
                );

                Info<< "    Number of faces shared with processor "
                    << ppp.neighbProcNo() << " = " << ppp.size() << endl;

                nProcPatches++;
                nProcFaces += ppp.size();
            }
            else
            {
                nBoundaryFaces += procMesh.boundaryMesh()[patchi].size();
            }
        }

        if (procMesh.nCells() && (nBoundaryFaces || nProcFaces))
        {
            Info<< "    Number of processor patches = " << nProcPatches << nl
                << "    Number of processor faces = " << nProcFaces << nl
                << "    Number of boundary faces = " << nBoundaryFaces << nl;
        }

        totProcFaces += nProcFaces;
        totProcPatches += nProcPatches;
        maxProcPatches = max(maxProcPatches, nProcPatches);
        maxProcFaces = max(maxProcFaces, nProcFaces);

        // create and write the addressing information
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procPointAddressing_[proci]
        );
        pointProcAddressing.write();

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procFaceAddressing_[proci]
        );
        faceProcAddressing.write();

        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procCellAddressing_[proci]
        );
        cellProcAddressing.write();

        // Write patch map for backwards compatibility.
        // (= identity map for original patches, -1 for processor patches)
        label nMeshPatches = curPatchSizes.size();
        labelList procBoundaryAddressing(identity(nMeshPatches));
        procBoundaryAddressing.setSize(nMeshPatches+nProcPatches, -1);

        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procBoundaryAddressing
        );
        boundaryProcAddressing.write();
    }

    scalar avgProcCells = scalar(nCells())/nProcs_;
    scalar avgProcPatches = scalar(totProcPatches)/nProcs_;
    scalar avgProcFaces = scalar(totProcFaces)/nProcs_;

    // In case of all faces on one processor. Just to avoid division by 0.
    if (totProcPatches == 0)
    {
        avgProcPatches = 1;
    }
    if (totProcFaces == 0)
    {
        avgProcFaces = 1;
    }

    Info<< nl
        << "Number of processor faces = " << totProcFaces/2 << nl
        << "Max number of cells = " << maxProcCells
        << " (" << 100.0*(maxProcCells-avgProcCells)/avgProcCells
        << "% above average " << avgProcCells << ")" << nl
        << "Max number of processor patches = " << maxProcPatches
        << " (" << 100.0*(maxProcPatches-avgProcPatches)/avgProcPatches
        << "% above average " << avgProcPatches << ")" << nl
        << "Max number of faces between processors = " << maxProcFaces
        << " (" << 100.0*(maxProcFaces-avgProcFaces)/avgProcFaces
        << "% above average " << avgProcFaces << ")" << nl
        << endl;

    return true;
}


// ************************************************************************* //
