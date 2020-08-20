    regionModels::surfaceFilmModels::SWIMIILayer* filmSWIM=dynamic_cast<regionModels::surfaceFilmModels::SWIMIILayer*> (&surfaceFilm);

    pointVectorField& pD=const_cast<pointVectorField&>(mesh.objectRegistry::lookupObject<pointVectorField>("pointDisplacement"));

    IOdictionary moveMeshesDict
    (
        IOobject
        (
            "moveMeshesDict",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );
    word mappedWallPatchName=word(moveMeshesDict.lookup("mappedWallPatchName")); //"region0_to_wallFilmRegion_wallFilmFaces";
    Info<<"Mapped patch: "<<mappedWallPatchName<<endl;
    label filmPatchID=mesh.boundaryMesh().findPatchID(mappedWallPatchName);
    vectorField& pDisp=refCast<vectorField>(pD.boundaryFieldRef()[filmPatchID]);

    fvMesh* rmesh = runTime.getObjectPtr<fvMesh>("wallFilmRegion");
    label intfilmPatchID=rmesh->boundaryMesh().findPatchID(mappedWallPatchName);

    word filmTopPatchName=word(moveMeshesDict.lookup("filmTopPatchName")); //"wallFilmFaces_top";
    Info<<"Film top patch: "<<filmTopPatchName<<endl;
    label inttopfilmPatchID=rmesh->boundaryMesh().findPatchID(filmTopPatchName);
    pointField newPs(rmesh->points().size());

/*
// cylinder radius
    scalar r=readScalar(moveMeshesDict.lookup("cylinderRadius")); //0.076; 
    Info<<"Cylinder radius: "<<r<<endl;
// film thickness (must be equal thickness in extrudeToRegionMeshDict)
    scalar filmThickness=readScalar(moveMeshesDict.lookup("filmThickness")); //0.005;
    Info<<"Film thickness: "<<filmThickness<<nl<<endl;
    scalar rtop=r-filmThickness;
*/    
    bool moveMeshes=readBool(moveMeshesDict.lookup("moveMeshes"));
    Info<<"Move meshes: "<<moveMeshes<<endl;
    
