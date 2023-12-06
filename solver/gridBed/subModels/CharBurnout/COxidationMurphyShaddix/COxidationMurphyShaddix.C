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

\*---------------------------------------------------------------------------*/

#include "COxidationMurphyShaddix.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class BedType>
Foam::label Foam::COxidationMurphyShaddix<BedType>::maxIters_ = 1000;

template<class BedType>
Foam::scalar Foam::COxidationMurphyShaddix<BedType>::tolerance_ = 1e-06;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::COxidationMurphyShaddix<BedType>::COxidationMurphyShaddix
(
    const dictionary& dict,
    BedType& owner
)
:
    SurfaceReactionModel<BedType>(dict, owner, typeName),
    D0_(readScalar(this->coeffDict().lookup("D0"))),
    rho0_(readScalar(this->coeffDict().lookup("rho0"))),
    T0_(readScalar(this->coeffDict().lookup("T0"))),
    Dn_(readScalar(this->coeffDict().lookup("Dn"))),
    A_(readScalar(this->coeffDict().lookup("A"))),
    E_(readScalar(this->coeffDict().lookup("E"))),
    n_(readScalar(this->coeffDict().lookup("n"))),
    WVol_(readScalar(this->coeffDict().lookup("WVol"))),
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
Foam::COxidationMurphyShaddix<BedType>::COxidationMurphyShaddix
(
    const COxidationMurphyShaddix<BedType>& srm
)
:
    SurfaceReactionModel<BedType>(srm),
    D0_(srm.D0_),
    rho0_(srm.rho0_),
    T0_(srm.T0_),
    Dn_(srm.Dn_),
    A_(srm.A_),
    E_(srm.E_),
    n_(srm.n_),
    WVol_(srm.WVol_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::COxidationMurphyShaddix<BedType>::~COxidationMurphyShaddix()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
Foam::scalar Foam::COxidationMurphyShaddix<BedType>::calculate
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
    // Surface combustion until combustible fraction is consumed
    if (mass < small)
    {
        return 0.0;
    }

    const SLGThermo& thermo = this->owner().thermo();

    // Cell carrier phase O2 species density [kg/m^3]
    const scalar rhoO2 = rhoc*thermo.carrier().Y(O2GlobalId_)[celli];

    if (rhoO2 < small)
    {
        return 0.0;
    }

    // Particle surface area [m^2]
    const scalar Ap = constant::mathematical::pi*sqr(d);

    // Calculate diffision constant at continuous phase temperature
    // and density [m^2/s]
    const scalar D = D0_*(rho0_/rhoc)*pow(Tc/T0_, Dn_);

    // Far field partial pressure O2 [Pa]
    const scalar ppO2 = rhoO2/WO2_*RR*Tc;

    // Total molar concentration of the carrier phase [kmol/m^3]
    const scalar C = pc/(RR*Tc);

    if (debug)
    {
        Pout<< "mass  = " << mass << nl
            << "Ap    = " << Ap << nl
            << "dt    = " << dt << nl
            << "C     = " << C << nl
            << endl;
    }

    // Molar reaction rate per unit surface area [kmol/m^2/s]
    scalar qCsOld = 0;
    scalar qCs = 1;

    const scalar qCsLim = mass/(WC_*Ap*dt);

    if (debug)
    {
        Pout<< "qCsLim = " << qCsLim << endl;
    }

    label iter = 0;
    while ((mag(qCs - qCsOld)/qCs > tolerance_) && (iter <= maxIters_))
    {
        qCsOld = qCs;
        const scalar PO2Surface = ppO2*exp(-(qCs + N)*d/(2*C*D));
        qCs = A_*exp(-E_/(RR*T))*pow(PO2Surface, n_);
        qCs = (100.0*qCs + iter*qCsOld)/(100.0 + iter);
        qCs = min(qCs, qCsLim);

        if (debug)
        {
            Pout<< "iter = " << iter
                << ", qCsOld = " << qCsOld
                << ", qCs = " << qCs
                << nl << endl;
        }

        iter++;
    }

    if (iter > maxIters_)
    {
        WarningInFunction
            << "iter limit reached (" << maxIters_ << ")" << nl << endl;
    }

    // Calculate the number of molar units reacted
    scalar dOmega = qCs*Ap*dt;

    // Add to carrier phase mass transfer
    dMassSRCarrier[O2GlobalId_] += -dOmega*WO2_;
    dMassSRCarrier[CO2GlobalId_] += dOmega*(WC_ + WO2_);

    // Add to particle mass transfer
    dmass += dOmega*WC_;

    const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // carrier sensible enthalpy exchange handled via change in mass

    // Heat of reaction [J]
    return dOmega*(WC_*HsC - (WC_ + WO2_)*HcCO2_);
}


// ************************************************************************* //
