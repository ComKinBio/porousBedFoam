/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    coalChemistryFoam

Description
    Transient solver for compressible, turbulent flow, with coal and limestone
    particle clouds, an energy source, and combustion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "coalChemistryTurbulenceModel.H"
#include "basicThermoCloud.H"
#include "coalCloud.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "fvOptions.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"

    //for DEM
    
    label DEMStep = 0;
    scalar DEMNextTime = 0.0;    
    
    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        //for DEM
        if (DEMSwitch == true)
        {
            if (coalParcels.DEMFlag == 1)
            {
                if (DEMStep>=DEMStepMax)
                {
                    coalParcels.DEMFlag = 0;
                    DEMStep = 0;
                    DEMNextTime = runTime.value() + DEMInter;
                }
                    
                DEMStep = DEMStep + 1;
            }
            else
            {
                if (runTime.value() >= DEMNextTime)
                {
                    coalParcels.DEMFlag = 1;
                }
            }
        }
        
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        rhoEffLagrangian = coalParcels.rhoEff() /*+ limestoneParcels.rhoEff()*/;
        pDyn = 0.5*rhoc*magSqr(Uc);

        coalParcels.evolve();

//         limestoneParcels.evolve();


        // Update continuous phase volume fraction field
        alphac = max(1.0 - coalParcels.theta() /*- limestoneParcels.theta()*/, alphacMin);
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
        //phic = linearInterpolate(Uc) & mesh.Sf();
        //rhocPhic = fvc::interpolate(rhoc)*phic;
        alphaRhoPhic = alphacf*rhocPhic;

        #include "rhocEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UcEqn.H"
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

        rhoc = thermo.rho();

        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
