/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "COxidationKDSteamAshMelt.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::COxidationKDSteamAshMelt<BedType>::
COxidationKDSteamAshMelt
(
    const dictionary& dict,
    BedType& owner
)
:
    SurfaceReactionModel<BedType>(dict, owner, typeName),
    COmega_(readScalar(this->coeffDict().lookup("COmiga"))),
    ep3_(readScalar(this->coeffDict().lookup("ep3"))),
    epLowT_(this->coeffDict().lookupOrDefault("epLowT", 0.9)),
    epHightT_(this->coeffDict().lookupOrDefault("epHightT", 0.01)),
    Tswitch_(this->coeffDict().lookupOrDefault("Tswitch", 1473.0)),
    delta_(this->coeffDict().lookupOrDefault("delta", 50.0)),
    C1_O2_(readScalar(this->coeffDict().lookup("C1_O2"))),
    C1_H2O_(readScalar(this->coeffDict().lookup("C1_H2O"))),
    C1_CO2_(readScalar(this->coeffDict().lookup("C1_CO2"))),
    C1_H2_(readScalar(this->coeffDict().lookup("C1_H2"))),
    C2_1_(readScalar(this->coeffDict().lookup("C2_1"))),
    C2_2_(readScalar(this->coeffDict().lookup("C2_2"))),
    C2_3_(readScalar(this->coeffDict().lookup("C2_3"))),
    C2_4_(readScalar(this->coeffDict().lookup("C2_4"))),
    E_1_(readScalar(this->coeffDict().lookup("E1"))),
    E_2_(readScalar(this->coeffDict().lookup("E2"))),
    E_3_(readScalar(this->coeffDict().lookup("E3"))),
    E_4_(readScalar(this->coeffDict().lookup("E4"))),
    CsLocalId_(-1),
    O2GlobalId_(owner.thermo().carrierId("O2")),
    CH4GlobalId_(owner.thermo().carrierId("CH4")),
    H2OGlobalId_(owner.thermo().carrierId("H2O")),
    CO2GlobalId_(owner.thermo().carrierId("CO2")),
    COGlobalId_(owner.thermo().carrierId("CO")),
    H2GlobalId_(owner.thermo().carrierId("H2")),
    WC_(0.0),
    WH2O_(0.0),
    WCO2_(0.0),
    WH2_(0.0),
    WO2_(0.0),
    HcCO2_(0.0),
    HcCO_(0.0)
{
    // Determine Cs ids
    CsLocalId_ = owner.thermo().solidId("C");

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().Wi(O2GlobalId_);
    WCO2_ = owner.thermo().carrier().Wi(CO2GlobalId_);
    WC_ = WCO2_ - WO2_;
    WH2O_ = owner.thermo().carrier().Wi(H2OGlobalId_);
    WH2_ = WH2O_ - WO2_/2.0;

    HcCO2_ = owner.thermo().carrier().Hc(CO2GlobalId_);
    HcCO_ = owner.thermo().carrier().Hc(COGlobalId_);
}


template<class BedType>
Foam::COxidationKDSteamAshMelt<BedType>::
COxidationKDSteamAshMelt
(
    const COxidationKDSteamAshMelt<BedType>& srm
)
:
    SurfaceReactionModel<BedType>(srm),
    //Sb_(srm.Sb_),
    COmega_(srm.COmega_),
    ep3_(srm.ep3_),
    epLowT_(srm.epLowT_),
    epHightT_(srm.epHightT_),
    Tswitch_(srm.Tswitch_),
    delta_(srm.delta_),
    C1_O2_(srm.C1_O2_),
    C1_H2O_(srm.C1_H2O_),
    C1_CO2_(srm.C1_CO2_),
    C1_H2_(srm.C1_H2_),
    C2_1_(srm.C2_1_),
    C2_2_(srm.C2_2_),
    C2_3_(srm.C2_3_),
    C2_4_(srm.C2_4_),
    E_1_(srm.E_1_),
    E_2_(srm.E_2_),
    E_3_(srm.E_3_),
    E_4_(srm.E_4_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CH4GlobalId_(srm.CH4GlobalId_),
    H2OGlobalId_(srm.H2OGlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    COGlobalId_(srm.COGlobalId_),
    H2GlobalId_(srm.H2GlobalId_),
    WC_(srm.WC_),
    WH2O_(srm.WH2O_),
    WCO2_(srm.WCO2_),
    WH2_(srm.WH2_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_),
    HcCO_(srm.HcCO_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::COxidationKDSteamAshMelt<BedType>::
~COxidationKDSteamAshMelt()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
Foam::scalar Foam::COxidationKDSteamAshMelt<BedType>::calculate
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
    const scalar massAsh,
    const scalar N,
    scalar& dmass,
    scalarField& dMassSRCarrier
) const
{
    // Surface combustion active combustible fraction is consumed
    if (mass < small)
    {
        return 0.0;
    }

    const SLGThermo& thermo = this->owner().thermo();
    
    ;
    const scalar YO2 = thermo.carrier().Y(O2GlobalId_)[celli];
    const scalar YCO2 = thermo.carrier().Y(CO2GlobalId_)[celli];
    const scalar YH2O = thermo.carrier().Y(H2OGlobalId_)[celli];
    const scalar YH2 = thermo.carrier().Y(H2GlobalId_)[celli];
    
    const scalar epsilon_ash = epLowT_ - (epLowT_ - epHightT_)/(1.0 + exp(-(T-Tswitch_)/delta_));

    // Local mole concentration of O2 in the carrier phase,            
    const scalar C_O2 = YO2*rhoc/WO2_;
    const scalar C_CO2 = YCO2*rhoc/WCO2_;
    const scalar C_H2O = YH2O*rhoc/WH2O_;
    const scalar C_H2 = YH2*rhoc/WH2_;

    // Diffusion rate coefficient, use O2 for all species
    const scalar D0 = C1_O2_*pow(0.5*(T + Tc), 2.0);
    const scalar D0ash = C1_O2_*pow(T, 2.0);

    // Schmidt number
    const scalar Sc = muc/(rhoc*D0);

    // Sherwood number
    const scalar Sh = this->owner().massTransfer().Sh(Rec, Sc);


    // Mass convection coefficient hma
    const scalar hma = D0*Sh/(d);
    
    // Assume ash intrinc density is 2000;
    const scalar dAsh = cbrt(6.0*massAsh/(constant::mathematical::pi*2000.0*(1-epsilon_ash)) + pow3(d));
    
    const scalar deltaD = max((dAsh - d)/2.0, small);

    // Mass diffusion coefficient through the ash layer hmi
    const scalar hmi = D0ash*pow(epsilon_ash, 2.0)/(deltaD);

    // General mass transfer rate for combustion
    const scalar beta_d = hmi*hma/(hmi+hma);

    // Kinetic rate beta_r
    const scalar beta_r_O2 = C2_1_*T*exp(-E_1_/T);
    const scalar beta_r_CO2 = C2_2_*T*exp(-E_2_/T);
    const scalar beta_r_H2O = C2_3_*T*exp(-E_3_/T);
    const scalar beta_r_H2 = C2_4_*T*exp(-E_4_/T);

    // Correlation for the formation of CO and CO2 omega_c
    const scalar omega_c = 2.0*(1.0+4.3*exp(-COmega_/T))/(2.0+4.3*exp(-COmega_/T));
    
    // Particle surface area for the char
    const scalar Ap = constant::mathematical::pi*sqr(d2nd);
    

    // Change in C mass [kg]
    scalar dmC_O2 = WC_*omega_c*C_O2*Ap*(beta_r_O2*beta_d/(beta_r_O2+beta_d))*dt;
    scalar dmC_CO2 = WC_*C_CO2*Ap*(beta_r_CO2*beta_d/(beta_r_CO2+beta_d))*dt;
    scalar dmC_H2O = WC_*C_H2O*Ap*(beta_r_H2O*beta_d/(beta_r_H2O+beta_d))*dt;
    scalar dmC_H2 = WC_*0.5*C_H2*Ap*(beta_r_H2*beta_d/(beta_r_H2+beta_d))*dt;

    // Limit mass transfer by availability of C
    scalar dmCTotal = dmC_O2 + dmC_CO2 + dmC_H2O + dmC_H2;
    dmC_O2 = min(mass*dmC_O2/dmCTotal, dmC_O2);
    dmC_CO2 = min(mass*dmC_CO2/dmCTotal, dmC_CO2);
    dmC_H2O = min(mass*dmC_H2O/dmCTotal, dmC_H2O);
    dmC_H2 = min(mass*dmC_H2/dmCTotal, dmC_H2);
    
    // Update local particle C mass
    dmass = dmC_O2 + dmC_CO2 + dmC_H2O + dmC_H2;

    // Molar consumption
    const scalar dOmega_O2 = dmC_O2/omega_c/WC_;
    const scalar dOmega_CO2 = dmC_CO2/WC_;
    const scalar dOmega_H2O = dmC_H2O/WC_;
    const scalar dOmega_H2 = dmC_H2/WC_*2.0;

    // Change in O2 mass [kg]
    const scalar dmO2 = dOmega_O2*WO2_;

    // Mass of newly created CO2 [kg]
    const scalar dmCO2 = (dOmega_O2*(2.0-omega_c)-dOmega_CO2)*(WC_ + WO2_);

    // Mass of newly created CO [kg]
    const scalar dmCO = (dOmega_O2*(2.0*(omega_c-1.0))+dOmega_CO2*2.0+dOmega_H2O)*(WC_ + WO2_/2.0);
    
    // Mass of newly created H2 [kg]
    const scalar dmH2 = (dOmega_H2O-dOmega_H2)*WH2_;
    
    // Mass of consumed H2O [kg]
    const scalar dmH2O = dOmega_H2O*WH2O_;
    
    // Mass of newly created CH4 [kg]
    const scalar dmCH4 = dOmega_H2*(WC_ + 2.0*WH2_)*0.5;

    
    // Update carrier O2 and CO2 mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] += dmCO2;
    dMassSRCarrier[COGlobalId_] += dmCO;
    dMassSRCarrier[H2OGlobalId_] -= dmH2O;
    dMassSRCarrier[CH4GlobalId_] += dmCH4;
    dMassSRCarrier[H2GlobalId_] += dmH2;
    
    const scalar LHV_index_CO2 = 0.0;               /* per definition */
    const scalar LHV_index_CO = 10.25e6;            /* J/kg D&D p. 26 */
    const scalar LHV_index_C = 32.79e6;         /* J/kg Alberto's code */
    const scalar LHV_index_CH4 = 50.0e6;         /* J/kg D&D p. 26 */
    const scalar LHV_index_H2 = 120.0e6;         /* J/kg D&D p. 26 */
//     const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // carrier sensible enthalpy exchange handled via change in mass
    // Heat of reaction [J]
    const scalar Heat_O2 = (dOmega_O2*(2.0-omega_c)*(WC_ + WO2_))*LHV_index_CO2\
                         + (dOmega_O2*(2.0*(omega_c-1.0))*(WC_ + WO2_/2.0))*LHV_index_CO - dmC_O2*LHV_index_C;
    const scalar Heat_CO2 = dOmega_CO2*(WC_ + WO2_/2.0)*2.0*LHV_index_CO-dOmega_CO2*LHV_index_C*WC_;
    const scalar Heat_H2O = dOmega_H2O*(LHV_index_H2*(WH2O_ - WO2_/2.0)+LHV_index_CO*(WC_ + WO2_/2.0)-LHV_index_C*WC_);
    const scalar Heat_H2 = 0.5*LHV_index_CH4*dOmega_H2*(WC_ + 2.0*(WH2O_ -WO2_/2.0))\
                          -dOmega_H2*LHV_index_H2*(WH2O_ - WO2_/2.0)-0.5*dOmega_H2*WC_*LHV_index_C;

    return -(Heat_O2 + Heat_CO2 + Heat_H2O + Heat_H2);
}


// ************************************************************************* //
