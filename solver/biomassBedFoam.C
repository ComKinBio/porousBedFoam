/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    biomassBedFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "bedTurbulenceModel.H"
#include "basicGravityBio2DBed.H"
#include "ParticleForce.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "fvOptions.H"
#include "SLGThermo.H"
#include "Random.H"
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
    
    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "\nStarting time loop\n" << endl;
    
    bool solverFirstIter = true;
 
    while (runTime.run())
    { 
    
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
        
        pDyn = 0.5*rhoc*magSqr(Uc);
       
        Info<< "Solve the bed model"
            << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
        bioBed.solveConversion();
        
        Info<< "Finish Solving the bed model"
            << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

            
        // update bed phase fraction
        alphac = max(1.0 - bioBed.beta(), alphacMin);
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
        alphaRhoPhic = alphacf*rhocPhic;
        
        if (bioBed.radiation() && 
            (solverFirstIter || (runTime.timeIndex() % radiationsolverFreq == 0)))
        {
            solverFirstIter = false;
            
            bedap = bioBed.ap();
            bedEp  = bioBed.Ep();
            bedsigmap = bioBed.sigmap();
        }
        
        
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
        
        Info<< "For one fluid timestep, ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
    }  
        
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
