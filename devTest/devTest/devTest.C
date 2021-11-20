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
    devTest

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "basic2DGridBed.H"
#include "basicGravity2DBed.H"
// #include "basicBio2DBed.H"
#include "ParticleForce.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "fvOptions.H"
#include "SLGThermo.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
        
    Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    SLGThermo slgThermo(mesh, thermo);

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();
        
    volVectorField Uc
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength/dimTime, vector::zero)
    );   
    
    
    
    surfaceScalarField phic
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Uc) & mesh.Sf()
    );

    volScalarField scalarFieldEmpty
    (
        IOobject
        (
            "alphac",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alphac",dimless, 0.0)
    );
    
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho",dimDensity, 1.1)
    );
    
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu",dimDensity, 0.00001)
    );
        
    
    const liquidMixtureProperties& liquids =  slgThermo.liquids();
    
    const scalar Dab = liquids.properties()[0].D(100000, 298);
    
    Info<<"Dab: "<<Dab<<nl<<endl;
    
    Info<<"slgThermo.liqidId(H2O): "<<slgThermo.liquidId("H2O")<<nl<<endl;
    
    Info<<"slgThermo.liquids().size(): "<<slgThermo.liquids().size()<<nl<<endl;
    
     Info<<"slgThermo.liquids().components(): "<<slgThermo.liquids().components()<<nl<<endl;
    
    Info<<"slgThermo.carrier().species().size(): "<<slgThermo.carrier().species().size()<<nl<<endl;
    
    Info<<"slgThermo.carrierId(H2O): "<<slgThermo.carrierId("H2O")<<nl<<endl;
    
    Info<<"slgThermo.carrierId(slgThermo.liquids().components()[0]): "<<slgThermo.carrierId(slgThermo.liquids().components()[0])<<nl<<endl;
    
    const label seed = 233;
    Random rndGen(12);
    
    Info<<" rndGen.sample01<scalar>(): "<< rndGen.sample01<scalar>()<<nl<<endl;
    Info<<" rndGen.sample01<scalar>(): "<< rndGen.sample01<scalar>()<<nl<<endl;
    Info<<" rndGen.sample01<scalar>(): "<< rndGen.sample01<scalar>()<<nl<<endl;
    Info<<" rndGen.sample01<scalar>(): "<< rndGen.sample01<scalar>()<<nl<<endl;
    Info<<" rndGen.sample01<scalar>(): "<< rndGen.sample01<scalar>()<<nl<<endl;
    
    
    
    Info<<" rndGen.sampleAB(1,23): "<< rndGen.sampleAB(1,23)<<nl<<endl;
    Info<<" rndGen.sampleAB(1,23): "<< rndGen.sampleAB(1,23)<<nl<<endl;
    Info<<" rndGen.sampleAB(1,23): "<< rndGen.sampleAB(1,23)<<nl<<endl;
    Info<<" rndGen.sampleAB(1,23): "<< rndGen.sampleAB(1,23)<<nl<<endl;
    Info<<" rndGen.sampleAB(1,23): "<< rndGen.sampleAB(1,23)<<nl<<endl;
    
    
    
   
    

//     basicBio2DBed testBed
//     (
//         "bed2d",
//         mesh,
//         rho,
//         Uc,
//         g,
//         slgThermo
//     );
    
//     basicGravity2DBed testBed
//     (
//         "bed2d",
//         mesh,
//         rho,
//         Uc,
//         g,
//         slgThermo
//     );
//     
//     const scalar dt = 0.0001;
//     
//     volScalarField& alphaBed = testBed.alpha();
//     
//     Info<<"alphaBed field: "<<alphaBed<<nl<<endl;
//     
//     Info<<"bed number field: "<<testBed.particleNumber()<<nl<<endl;
//     
//     testBed.gravityCollapse();
//     
//     Info<<"gravityCollapse: "<<nl<<endl;
//     
//     Info<<"alphaBed field: "<<alphaBed<<nl<<endl;
//     
//     Info<<"bed number field: "<<testBed.particleNumber()<<nl<<endl;
//     
//     testBed.solve(dt);
//     
//     fvVectorMatrix SU = testBed.SU(Uc,dt);
//     
// 
//     Info<< "End\n" << endl;
//     Info<< "SU" << SU << nl << endl;

    return 0;
}


// ************************************************************************* //
