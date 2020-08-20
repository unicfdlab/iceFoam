/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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
    reactingParcelSWIMFoam

Group
    grpLagrangianSolvers

Description
    Transient solver for compressible, turbulent flow with a reacting,
    multiphase particle cloud, and surface film modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "motionSolver.H"
#include "displacementLaplacianFvMotionSolver.H"
#include "turbulentFluidThermoModel.H"
#include "basicReactingMultiphaseCloud.H"
#include "surfaceFilmModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "SWIMIILayer.H"
#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"
#include "volPointInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])

{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUf.H"
    #include "createRegionControls.H"
    #include "initContinuityErrs.H"
    #include "initMoveMeshes.H"


    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    volScalarField olddeltaIce=filmSWIM->deltaIce();

    while (runTime.run())
    {
        #include "readDyMControls.H"


        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setMultiRegionDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();
	TParcels=dimensionedScalar("", dimTemperature, parcels.Tmin());
	olddeltaIce=filmSWIM->deltaIce();
        surfaceFilm.evolve();
        
	if (moveMeshes) {
            #include "moveMeshes.H"
        }

        if (solvePrimaryRegion)
        {
            if (pimple.nCorrPIMPLE() <= 1)
            {
                #include "rhoEqn.H"
            }

            // --- PIMPLE loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YEqn.H"
                #include "EEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
