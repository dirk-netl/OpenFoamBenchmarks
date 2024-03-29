    Info<< "Create dynamic mesh for time = "
        << runTime.timeName() << nl << endl;

    autoPtr<dynamicFvMesh> meshPtr
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    dynamicFvMesh& mesh = meshPtr();

    // Calculate initial mesh-to-mesh mapping. Note that this should be
    // done under the hood, e.g. as a MeshObject
    mesh.update();
