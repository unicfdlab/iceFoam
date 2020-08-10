# iceFoam

UniCFD lab solver for ice accretion simulation.

Current status of folders and files.

├── iceFoam - solver
│       ├── createClouds.H
│       ├── createFieldRefs.H
│       ├── createFields.H
│       ├── createRegionControls.H
│       ├── createSurfaceFilmModel.H
│       ├── EEqn.H
│       ├── Make
│       │   ├── files
│       │   └── options
│       ├── pEqn.H
│       ├── iceFoam.C
│       ├── rhoEqn.H
│       ├── setMultiRegionDeltaT.H
│       ├── setRDeltaT.H
│       ├── UEqn.H
│       └── YEqn.H
│
├── subModels
│   ├── surfaceFilmModelsSWIM - SWIM model
│   │   ├── Make
│   │   │   ├── files
│   │   │   └── options
│   │   ├── README.MD
│   │   └── SWIMIILayer
│   │       ├── SWIMIILayer.C
│   │       ├── SWIMIILayer.H
│   │       └── SWIMIILayerI.H
│   └── particlesLib - lib of extra particle features (in developing)
│
├── tutorials -- cases
│   ├── cylinder2D -- cylinder
│   └── NACA0012   -- NACA airfoil
│
└── README.MD
