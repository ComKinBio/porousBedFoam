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

\*---------------------------------------------------------------------------*/

#include "COxidationKineticDiffusionLimitedRate.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::COxidationKineticDiffusionLimitedRate<BedType>::
COxidationKineticDiffusionLimitedRate
(
    const dictionary& dict,
    BedType& owner
)
:
    SurfaceReactionModel<BedType>(dict, owner, typeName),
    Sb_(readScalar(this->coeffDict().lookup("Sb"))),
    C1_(readScalar(this->coeffDict().lookup("C1"))),
    C2_(readScalar(this->coeffDict().lookup("C2"))),
    E_(readScalar(this->coeffDict().lookup("E"))),
    CsLocalId_(-1),
    O2GlobalId_(owner.thermo().carrierId("O2")),
    CO2GlobalId_(owner.thermo().carrierId("CO2")),
    WC_(0.0),
    WO2_(0.0),
    HcCO2_(0.0)
{
    // Determine Cs ids
    CsLocalId_ = owner.thermo().solidId("C");

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().Wi(O2GlobalId_);
    const scalar WCO2 = owner.thermo().carrier().Wi(CO2GlobalId_);
    WC_ = WCO2 - WO2_;

    HcCO2_ = owner.thermo().carrier().Hc(CO2GlobalId_);
}


template<class BedType>
Foam::COxidationKineticDiffusionLimitedRate<BedType>::
COxidationKineticDiffusionLimitedRate
(
    const COxidationKineticDiffusionLimitedRate<BedType>& srm
)
:
    SurfaceReactionModel<BedType>(srm),
    Sb_(srm.Sb_),
    C1_(srm.C1_),
    C2_(srm.C2_),
    E_(srm.E_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::COxidationKineticDiffusionLimitedRate<BedType>::
~COxidationKineticDiffusionLimitedRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
Foam::scalar Foam::COxidationKineticDiffusionLimitedRate<BedType>::calculate
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

    // Local mass fraction of O2 in the carrier phase
    const scalar YO2 = thermo.carrier().Y(O2GlobalId_)[celli];

    // Diffusion rate coefficient
    const scalar D0 = C1_/d*pow(0.5*(T + Tc), 0.75);

    // Kinetic rate
    const scalar Rk = C2_*exp(-E_/(RR*Tc));

    // Particle surface area
    const scalar Ap = constant::mathematical::pi*sqr(d);

    // Change in C mass [kg]
    scalar dmC = Ap*rhoc*RR*Tc*YO2/WO2_*D0*Rk/(D0 + Rk)*dt;

    // Limit mass transfer by availability of C
    dmC = min(mass, dmC);

    // Molar consumption
    const scalar dOmega = dmC/WC_;

    // Change in O2 mass [kg]
    const scalar dmO2 = dOmega*Sb_*WO2_;

    // Mass of newly created CO2 [kg]
    const scalar dmCO2 = dOmega*(WC_ + Sb_*WO2_);

    // Update local particle C mass
    dmass += dOmega*WC_;

    // Update carrier O2 and CO2 mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] += dmCO2;

    const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // carrier sensible enthalpy exchange handled via change in mass

    // Heat of reaction [J]
    return dmC*HsC - dmCO2*HcCO2_;
}


// ************************************************************************* //
