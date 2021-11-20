/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::BioBed

Description
    BioBed main functions, belongs to BioBed.C

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class BedType>
void Foam::BioBed<BedType>::solveConversion(const scalar dt)
{

    preSolve();
    
    tmp<volScalarField> Gc
    (
        volScalarField::New
        (
            "G_coppy",
            this->mesh_,
            dimensionedScalar(dimPower/sqr(dimLength), 0)
        )
    ); 
    
    if (this->radiation())
    {
        Gc.ref() = this->mesh_.objectRegistry::template
                    lookupObject<volScalarField>("G");
    }
    
    const labelList bedList = this->bedIDList();
    //- In submodels that not really used in the calculation
    scalar defaultScalar = 1;
    const scalar defaultConstScalar = 1;
    
     // Define local properties at beginning of timestep
    volScalarField& npWet = particleNumbere(wet_);
    volScalarField& npDry = particleNumbere(dry_);
    volScalarField& npChar = particleNumbere(char_);
    volScalarField& npAsh = particleNumbere(ash_);
    
    volScalarField& dWet = dpe(wet_);
    volScalarField& dDry = dpe(dry_);
    volScalarField& dChar = dpe(char_);
    volScalarField& dAsh = dpe(ash_);
    
    volScalarField& d2ndWet = dp2nde(wet_);
    volScalarField& d2ndDry = dp2nde(dry_);
    volScalarField& d2ndChar = dp2nde(char_);
    volScalarField& d2ndAsh = dp2nde(ash_);
    
    volScalarField& T0Wet = bedTe(wet_);
    volScalarField& T0Dry = bedTe(dry_);
    volScalarField& T0Char = bedTe(char_);
    volScalarField& T0Ash = bedTe(ash_);
    
    volScalarField& mass0Wet = masse(wet_);
    volScalarField& mass0Dry = masse(dry_);
    volScalarField& mass0Char = masse(char_);
    volScalarField& mass0Ash = masse(ash_);
    
    volScalarField& cpWet = bedCpe(wet_);
    volScalarField& cpDry = bedCpe(dry_);
    volScalarField& cpChar = bedCpe(char_);
    volScalarField& cpAsh = bedCpe(ash_);
    
    volScalarField& w = w_percent();
    volScalarField& gamma = gamma_percent();
    volScalarField& eta = eta_percent();
    
    const scalar epsilon0Wet = epsilon0e(wet_);
    const scalar epsilon0Dry = epsilon0e(dry_);
    const scalar epsilon0Char = epsilon0e(char_);
    const scalar epsilon0Ash = epsilon0e(ash_);
    
    const scalar kappaWet = kp0e(wet_);
    const scalar kappaDry = kp0e(dry_);
    const scalar kappaChar = kp0e(char_);
    const scalar kappaAsh = kp0e(ash_);

    //- bed loop
    forAll(bedList, i)
    {
        const label celli = bedList[i];
        const scalar pc = this->pc()[i];
        const scalar Tc = this->Tc()[i];
        const scalar Uc = this->U()[i];
        scalar G = 0.0;
        
        if (this->radiation())
        {
            G = Gc()[celli];
        }  
        
        // Sources
        //~~~~~~~~

        // Explicit momentum source for particle
        vector Su = Zero;

        // Linearised momentum source coefficient
        scalar Spu = 0.0;

        // Momentum transfer from the particle to the carrier phase
        vector dUTrans = Zero;

        // Explicit enthalpy source for particle
        scalar She = 0.0;

        // Linearised enthalpy source coefficient
        scalar Sphe = 0.0;

        // Sensible enthalpy transfer from the particle to the carrier phase
        scalar dhsTrans = 0.0;
        
        // Molar flux of species emitted from the particle (kmol/m^2/s)
        scalar Ne = 0.0;

        // Sum Ni*Cpi*Wi of emission species
        scalar NCpW = 0.0;
        
        // Surface concentrations of emitted species
        scalarField Cs(this->thermo().carrier().species().size(), 0.0);
        
        // Suface thermo veriables
        scalar Ts, rhos, mus, Prs, kappas;

        // 1. Compute models that contribute to mass transfer - U, T held constant
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            // Phase change in liquid phase
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            //TODO if water is the only liquids, code can be rewriten for water only
            // Mass transfer due to phase change
            scalarField dMassDrying(this->thermo().liquids().size(), 0.0);
            
            scalarField Yliquids(1,1);//Only one species, its mass fraction is 1
            
            const scalar npWeti = npWet[celli];
            if (npWeti >= 1)
            {
                const scalar dWeti = dWet[celli];
                const scalar d2ndWeti = d2ndWet[celli];
                const scalar T0Weti = T0Wet[celli];
                const scalar mass0Wateri = mass0Wet[celli]*w[celli];

                // Calc surface values for wet particle
                this->calcSurfaceValues(celli, T0Weti, Ts, rhos, mus, Prs, kappas);
                scalar Res = this->Re(rhos, Uc, d2ndWeti, mus);
                
                scalar shWet = 0.0;
                scalar sphWet = 0.0;
                scalar dhsTransWet = 0.0;

                // Calc mass and enthalpy transfer due to phase change
                calcDrying
                (
                    dt,
                    celli,
                    Res,
                    Prs,
                    defaultConstScalar,
                    defaultConstScalar,
                    mus/rhos,
                    d2ndWeti,
                    T0Weti,
                    Ts,
                    Tc,
                    pc,
                    mass0Wateri,
                    npWeti,
                    Yliquids,
                    dMassDrying,
                    shWet,
                    Ne,
                    NCpW,
                    Cs
                );
                
                // dMassDrying is the total mass of vaparized water
                // Update wet bed fields
                if (mag(dMassDrying - mass0Wateri)<small)
                {
                    T0Dry[celli] = (T0Dry[celli]*cpDry[celli]*mass0Dry[celli]*npDry[celli]
                                  + T0Wet[celli]*cpWet[celli]*mass0Wet[celli]*npWet[celli])\
                                  /(cpDry[celli]*mass0Dry[celli]*npDry[celli]
                                  + cpWet[celli]*mass0Wet[celli]*npWet[celli])
                    
                    dWet[celli] = dpAfterShrink(dWeti, mass0Wet[celli], mass0Wateri, shrinkageFactorAlpha());
                    d2ndWet[celli] = dWet[celli]/sqrt(this->sphericity());

                    dDry[celli] = cbrt((pow3(dDry[celli])*npDry[celli] + pow3(dWet[celli])*npWeti)\
                                 /(npDry[celli] + npWet[celli]));
                    d2ndDry[celli] = sqrt((pow3(d2ndDry[celli])*npDry[celli] + pow3(d2ndWet[celli])*npWeti)\
                                    /(npDry[celli] + npWet[celli]));
                    
                    gamma[celli] = mass0Dry[celli]*gamma[celli]*npDry[celli]\
                                  /(mass0Dry[celli]*(1 - gamma[celli])*npDry[celli]\
                                  + mass0Wet[celli](1- w[celli])*npWeti);
                    mass0Dry[celli] += mass0Wet[celli](1- w[celli]);
                    
                    w[celli] = 0;
                    mass0Wet[celli] = 0;
                    
                    npDry[celli] += npWet[celli];
                    npWet[celli] = 0;
                }
                else
                {
                    dWet[celli] = dpAfterShrink(dWeti, mass0Wet[celli], dMassDrying, shrinkageFactorAlpha());
                    d2ndWet[celli] = dWet[celli]/sqrt(this->sphericity());
                    mass0Wet[celli] -= dMassDrying;
                    w[celli] = (mass0Wateri - dMassDrying)/mass0Wet[celli];
                }
                
                // Calculate new particle temperature
                // Correct surface values due to emitted species
//                 this->correctSurfaceValues(Ts, Cs, rhos, mus, Prs, kappas);
                Res = this->Re(rhos, Uc, d2ndWet[celli], mus);
                
                // Calculate new particle temperature
                const scalar cpWater = cpWater[celli];
                
                T0Wet[celli] =
                    calcHeatTransfer
                    (
                        dt,
                        dWet[celli],
                        d2ndWet[celli],
                        mass0Wet[celli],
                        T0Wet[celli],
                        cpWater,
                        Tc,
                        Res,
                        Prs,
                        kappaWet,
                        NCpW,
                        shWet,
                        G,
                        epsilon0Wet,
                        dhsTransWet,
                        sphWet
                    );
                
               Sphe += sphWet*npWet[celli];
               dhsTrans += dhsTransWet*npWet[celli];
                        
            }
            
            

            // Devolatilisation
            // ~~~~~~~~~~~~~~~~
            
            // Mass transfer due to devolatilisation
            scalarField dMassDV(this->thermo().carrier().species().size(), 0.0);
            scalar dMass = 0;
            scalar dMassChar = 0;
                
            const scalar npDryi = npDry[celli];
            if (npDryi >= 1)
            {
                const scalar dDryi = dDry[celli];
                const scalar d2ndDryi = d2ndDry[celli];
                const scalar T0Dryi = T0Dry[celli];
                const scalar mass0YieldChari = mass0Dry[celli]*gamma[celli];
                const scalar mass0ActiveDryi = mass0Dry[celli]*(1-gamma[celli]);
            
                
                // Calc surface values for wet particle
                this->calcSurfaceValues(celli, T0Dryi, Ts, rhos, mus, Prs, kappas);
                
                scalar shDry = 0.0;
                scalar sphDry = 0.0;
                scalar dhsTransDry = 0.0;
                Ne = 0.0;
                NCpW = 0.0;

                // Calc mass and enthalpy transfer due to devolatilisation
                calcDevo
                (
                    dt,
                    defaultConstScalar,
                    Ts,
                    d2ndDryi,
                    T0Dryi,
                    mass0ActiveDryi,
                    mass0ActiveDryi,
                    npDryi,
                    rhos,
                    pc,
                    Tc,
                    true,
                    dMassDV,
                    dMass,
                    dMassChar,
                    shDry,
                    Ne,
                    NCpW,
                    Cs
                );
                
                // dMass is the total mass of devo wood
                if (mag(dMass - mass0ActiveDryi)<small)
                {
                    T0Char[celli] = (T0Char[celli]*cpChar[celli]*mass0Char[celli]*npChar[celli]
                                  + T0Dry[celli]*cpDry[celli]*mass0Dry[celli]*npDry[celli])\
                                  /(cpChar[celli]*mass0Char[celli]*npChar[celli]
                                  + cpDry[celli]*mass0Dry[celli]*npDry[celli])
                    
                    dDry[celli] = dpAfterShrink(dDryi, mass0Dry[celli], dMass, shrinkageFactorBeta());
                    d2ndDry[celli] = dDry[celli]/sqrt(this->sphericity());

                    dChar[celli] = cbrt((pow3(dChar[celli])*npChar[celli] + pow3(dDry[celli])*npDryi)\
                                 /(npChar[celli] + npDry[celli]));
                    d2ndChar[celli] = sqrt((pow3(d2ndChar[celli])*npChar[celli] + pow3(d2ndDry[celli])*npDryi)\
                                    /(npChar[celli] + npDry[celli]));
                    
                    eta[celli] = mass0Char[celli]*eta[celli]*npChar[celli]\
                                  /(mass0Char[celli]*(1 - eta[celli])*npChar[celli]\
                                  + mass0Dry[celli](1- gamma[celli])*npDryi);
                    mass0Char[celli] += mass0Dry[celli](1- gamma[celli]);
                    
                    gamma[celli] = 0;
                    mass0Dry[celli] = 0;
                    
                    npChar[celli] += npDry[celli];
                    npDry[celli] = 0;
                }
                else
                {
                    dDry[celli] = dpAfterShrink(dDryi, mass0Dry[celli], dMass, shrinkageFactorBeta());
                    d2ndDry[celli] = dDry[celli]/sqrt(this->sphericity());
                    gamma[celli] = mass0YieldChari/(mass0Dry[celli]- dMass);
                    mass0Dry[celli] -= dMass;
                }
                
                // Calculate new particle temperature
                // Correct surface values due to emitted species
//                 this->correctSurfaceValues(Ts, Cs, rhos, mus, Prs, kappas);
                Res = this->Re(rhos, Uc, d2ndDry[celli], mus);
                
                // Calculate new particle temperature
                const scalar cpDry = cpDry[celli];

                T0Dry[celli] =
                    calcHeatTransfer
                    (
                        dt,
                        dDry[celli],
                        d2ndDry[celli],
                        mass0Dry[celli],
                        T0Dry[celli],
                        cpDry,
                        Tc,
                        Res,
                        Prs,
                        kappaDry,
                        NCpW,
                        shDry,
                        G,
                        epsilon0Dry,
                        dhsTransDry,
                        sphDry
                    );
                
               Sphe += sphDry*npDry[celli];
               dhsTrans += dhsTransDry*npDry[celli];
            
            }
            
            
            // Surface reactions
            // ~~~~~~~~~~~~~~~~
            
            // Change in carrier phase composition due to surface reactions
            scalarField dMassSRCarrier(this->thermo().carrier().species().size(), 0.0);
            scalar dmassCharBurned = 0;
            
            const scalar npChari = npChar[celli];
            if (npChari >= 1)
            {
                const scalar dChari = dChar[celli];
                const scalar d2ndChari = d2ndChar[celli];
                const scalar T0Chari = T0Char[celli];
                const scalar mass0AshInChari = mass0Char[celli]*eta[celli];
                const scalar mass0ActiveChari = mass0Char[celli]*(1-eta[celli]);
            
                
                // Calc surface values for wet particle
                this->calcSurfaceValues(celli, T0Chari, Ts, rhos, mus, Prs, kappas);
                scalar Res = this->Re(rhos, Uc, d2ndChari, mus);
        //         scalar e_bed = this->alpha()[celli]
                
                scalar shChar = 0.0;
                scalar sphChar = 0.0;
                scalar dhsTransChar = 0.0;
                Ne = 0.0;
                NCpW = 0.0;

                
                // Calc mass and enthalpy transfer due to surface reactions
                calcCharBurnOut
                (
                    dt,
                    celli,
                    dChari,
                    d2ndChari,
                    T0Chari,
                    Tc,
                    pc,
                    rhos,
                    mus,
                    Res,
                    0.4, //e_bed, hardcoded now
                    mass0ActiveChari,
                    npChari,
                    Ne,
                    dmassCharBurned,
                    dMassSRCarrier,
                    shChar,
                    dhsTransChar
                );
                
                // dmassCharBurned is the total mass of char burnout
                if (mag(dmassCharBurned - mass0ActiveChari)<small)
                {
                    T0Ash[celli] = (T0Ash[celli]*cpAsh[celli]*mass0Ash[celli]*npAsh[celli]
                                  + T0Char[celli]*cpChar[celli]*mass0Char[celli]*npChar[celli])\
                                  /(cpAsh[celli]*mass0Ash[celli]*npAsh[celli]
                                  + cpChar[celli]*mass0Char[celli]*npChar[celli])
                    
                    mass0Ash[celli] += mass0Char[celli](1- gamma[celli]);
                    
                    eta[celli] = 0;
                    mass0Char[celli] = 0;
                    
                    npAsh[celli] += npChar[celli];
                    npChar[celli] = 0;
                }
                else
                {
                    dChar[celli] = dpAfterShrink(dChari, mass0Char[celli], dmassCharBurned, shrinkageFactorBeta());
                    d2ndChar[celli] = dChar[celli]/sqrt(this->sphericity());
                    eta[celli] = mass0AshInChari/(mass0Char[celli]- dmassCharBurned);
                    mass0Char[celli] -= dmassCharBurned;
                }
                
                // Calculate new particle temperature
                // Correct surface values due to emitted species
//                 this->correctSurfaceValues(Ts, Cs, rhos, mus, Prs, kappas);
                Res = this->Re(rhos, Uc, d2ndChar[celli], mus);
                
                // Calculate new particle temperature
                const scalar cpChar = cpChar[celli];
                
                T0Char[celli] =
                    calcHeatTransfer
                    (
                        dt,
                        dChar[celli],
                        d2ndChar[celli],
                        mass0Char[celli],
                        T0Char[celli],
                        cpChar,
                        Tc,
                        Res,
                        Prs,
                        kappaChar,
                        NCpW,
                        shChar,
                        G,
                        epsilon0Char,
                        dhsTransChar,
                        sphChar
                    );
                
               Sphe += sphChar*npChar[celli];
               dhsTrans += dhsTransChar*npChar[celli];
                
            }
            
            const scalar npAshi = npAsh[celli];
            if (npChari >= 1)
            {
                const scalar dAshi = dAsh[celli];
                const scalar d2ndAshi = d2ndAsh[celli];
                const scalar T0Ashi = T0Ash[celli];
                
                scalar sphAsh = 0.0;
                scalar dhsTransAsh = 0.0;
                
                Res = this->Re(rhos, Uc, d2ndAshi[celli], mus);
                
                // Calculate new particle temperature
                const scalar cpAsh = cpAsh[celli];
                
                T0Ash[celli] =
                    calcHeatTransfer
                    (
                        dt,
                        dAshi,
                        d2ndAshi,
                        mass0Ash[celli],
                        T0Ashi,
                        cpAsh,
                        Tc,
                        Res,
                        Prs,
                        kappaAsh,
                        0.0,
                        0.0,
                        G,
                        epsilon0Ash,
                        dhsTransAsh,
                        sphAsh
                    );
                
               Sphe += sphAsh*npAshi;
               dhsTrans += dhsTransAsh*npAshi;
            }

        // 2. Update particle average properites and calc momentum transfers
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Update particle average properites
        // ~~~~~~~~~~~~~
        
        
        // Correct surface values due to emitted species
//         this->correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);
        Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);
        
        
        // Motion
        // ~~~~~~

        // Calculate new particle velocity
        this->U_ =
            this->calcVelocity(cloud, td, dt, Res, mus, mass1, Su, dUTrans, Spu);


        // 3. Accumulate carrier phase source terms
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(dMassDrying, i)
        {
            const label cid = this->thermo().carrierId(this->thermo().liquids().components()[i]);
            scalar dm = npWeti*dMassDrying[i];
            scalar hs = composition.carrier().Hs(cid, pc, T0Weti);
            cloud.rhoTrans(cid)[celli] += dm;
            cloud.hsTrans()[celli] += dm*hs;
        }
        
        forAll(dMassDV, i)
        {
            scalar dm = npDryi*dMassDV[i];
            scalar hs = composition.carrier().Hs(i, pc, T0Dryi);
            cloud.rhoTrans(i)[this->cell()] += dm;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }
        
        forAll(dMassSRCarrier, i)
        {
            scalar dm = npChari*dMassSRCarrier[i];
            scalar hs = composition.carrier().Hs(i, pc, T0Chari);
            cloud.rhoTrans(i)[this->cell()] += dm;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }

//         // Update momentum transfer
//         cloud.UTrans()[this->cell()] += np0*dUTrans;
//         cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTrans()[celli] += dhsTrans;
        cloud.hsCoeff()[celli] += Sphe;

        // Update radiation fields
        if (this->radiation())
        {
            const scalar apWet = this->areaP(d2ndWet[celli]);
            const scalar T4Wet = pow4(T0Wet[celli]);
            const scalar apDry = this->areaP(d2ndDry[celli]);
            const scalar T4Dry = pow4(T0Dry[celli]);
            const scalar apChar = this->areaP(d2ndChar[celli]);
            const scalar T4Char = pow4(T0Char[celli]);
            const scalar apAsh = this->areaP(d2ndAsh[celli]);
            const scalar T4Ash = pow4(T0Ash[celli]);
            
            cloud.radAreaP()[celli] += 
                dt*(npWet[celli]*apWet + npDry[celli]*apDry + npChar[celli]*apChar + npAsh[celli]*apAsh);
            cloud.radT4()[celli] += 
                dt*(npWet[celli]*T4Wet + npDry[celli]*T4Dry + npChar[celli]*T4Char + npAsh[celli]*T4Ash);
            cloud.radAreaPT4()[celli] += 
                dt*(npWet[celli]*apWet*T4Wet + npDry[celli]*apDry*T4Dry + npChar[celli]*apChar*T4Char + npAsh[celli]*apAsh*T4Ash);
        }
        
    }//end bedIDloop
}


template<class BedType>
void Foam::BioBed<BedType>::calcDrying
(
    const scalar dt,
    const label celli,
    const scalar Re,
    const scalar Pr,
    const scalar Nu,
    const scalar Sh,
    const scalar nus,
    const scalar d,
    const scalar T,
    const scalar Ts,
    const scalar Tc,
    const scalar pc,
    const scalar mass,
    const scalar np,
    const scalarField& Y,
    scalarField& dMassPC,
    scalar& She,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
)
{
    scalarField X(this->thermo().liquids().X(Y));

    scalar Tvap = drying().Tvap(X);

    if (T < Tvap)
    {
        return;
    }

    const scalar TMax = drying().TMax(pc, X);
    const scalar Tdash = min(T, TMax);
    const scalar Tsdash = min(Ts, TMax);

    // Calculate mass transfer due to phase change
    drying().calculate
    (
        dt,
        celli,
        Re,//
        Nu,//
        Sh,
        d,
        nus,//
        Tdash,
        Tsdash,
        pc,
        Tc,
        X,
        dMassPC
    );

    // Limit phase change mass by availability of each specie
    dMassPC = min(mass, dMassPC);

    const scalar dMassTot = sum(dMassPC);

    // Add to cumulative phase change mass
    drying().addToPhaseChangeMass(np*dMassTot);
    
    
    forAll(dMassPC, i)
    {
        const label cid = this->thermo().carrierId(this->thermo().liquids().components()[i]);

        const scalar dh = drying().dh(cid, i, pc, Tdash);
        She -= dMassPC[i]*dh/dt;
    }


    // Update molar emissions
    if (this->heatTransfer().BirdCorrection())
    {
        // Average molecular weight of carrier mix - assumes perfect gas
        const scalar Wc = this->rho()[celli]*RR*Tc/pc;

        forAll(dMassPC, i)
        {
            const label cid = this->thermo().carrierId(this->thermo().liquids().components()[i]);
            
            const scalar Cp = this->thermo().carrier().Cp(cid, pc, Tsdash);
            const scalar W = this->thermo().carrier().Wi(cid);
            const scalar Ni = dMassPC[i]/(this->areaS(d)*dt*W);

            const scalar Dab =
                this->thermo().liquids().properties()[i].D(pc, Tsdash, Wc);

            // Molar flux of species coming from the particle (kmol/m^2/s)
            N += Ni;

            // Sum of Ni*Cpi*Wi of emission species
            NCpW += Ni*Cp*W;

            // Concentrations of emission species
            Cs[cid] += Ni*d/(2.0*Dab);
        }
    }
}


template<class BedType>
void Foam::BioBed<BedType>::calcDevo
(
    const scalar dt,
    const scalar age,
    const scalar Ts,
    const scalar d2nd,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const scalar np,
    const scalar rhoc,
    const scalar pc,
    const scalar Tc,
    label& canCombust,
    scalarField& dMassDV,
    scalar& dMass,
    scalar& dMassChar,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
)
{
   
    // Initialise demand-driven constants
    (void)TDevol();
    (void)LDevol();

    // Check that the parcel temperature is within necessary limits for
    // devolatilisation to occur
    if (T < TDevol())
    {
        return;
    }

    // Total mass of volatiles evolved
    devolatilisation().calculate
    (
        dt,
        age,//age no use
        mass0,//no use now
        mass,
        T,
        canCombust,//no use now
        dMass,
        dMassDV,
        dMassChar
    );

    scalar dMassTot = sum(dMassDV);

    devolatilisation().addToDevolatilisationMass
    (
        np*dMassTot
    );

    Sh -= dMassTot*LDevol()/dt;

    // Update molar emissions
    if (this->heatTransfer().BirdCorrection())
    {
        // Molar average molecular weight of carrier mix
        const scalar Wc = max(small, rhoc*RR*Tc/pc);

        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));

        forAll(dMassDV, i)
        {
            // i is already global gas species id
            const scalar Cp = this->thermo.carrier().Cp(i, pc, Ts);
            const scalar W = this->thermo.carrier().Wi(i);
            const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);

            // Dab calc'd using API vapour mass diffusivity function
            const scalar Dab =
                3.6059e-3*(pow(1.8*Ts, 1.75))
               *sqrt(1.0/W + 1.0/Wc)
               /(pc*beta);

            N += Ni;
            NCpW += Ni*Cp*W;
            Cs[i] += Ni*d/(2.0*Dab);
        }
    }
}


template<class BedType>
void Foam::BioBed<BedType>::calcCharBurnOut
(
    const scalar dt,
    const label celli,
    const scalar d,
    const scalar d2nd,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar muc,
    const scalar Rec,
    const scalar e_bed,
    const scalar mass,
    const scalar np,
    const scalar N,
    scalar& dmass,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans
)
{

    // Initialise demand-driven constants
    (void)hRetentionCoeff();
    (void)this->TMax();

    // Update surface reactions
    const scalar hReaction = surfaceReaction().calculate
    (
        dt,
        celli,
        d,
        d2nd,
        T,
        Tc,
        pc,
        rhoc,
        muc,
        Rec,
        e_bed,
        mass,
        N,
        dmass,
        dMassSRCarrier
    );

    surfaceReaction().addToSurfaceReactionMass(np*dmass);

    const scalar xsi = min(T/this->TMax(), 1.0);
    const scalar coeff =
        (1.0 - xsi*xsi)*hRetentionCoeff();

    Sh += coeff*hReaction/dt;

    dhsTrans += (1.0 - coeff)*hReaction;
}

// ************************************************************************* //
