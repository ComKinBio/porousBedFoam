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

#include "ThermoBed.H"
#include "meshTools.H"

#include "HeatTransferModel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class BedType>
void Foam::ThermoBed<BedType>::setModels()
{
    heatTransferModel_.reset
    (
        HeatTransferModel<ThermoBed<BedType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
   
    this->subModelProperties().lookup("radiation") >> radiation_;

    if (radiation_)
    {
        radAreaP_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->bedName() + ":radAreaP",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimArea, 0)
            )
        );

        radT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->bedName() + ":radT4",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(pow4(dimTemperature), 0)
            )
        );

        radAreaPT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    this->bedName() + ":radAreaPT4",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(sqr(dimLength)*pow4(dimTemperature), 0)
            )
        );
    }
    
    dictionary integrationSchemeDict_ = this->subModelProperties_.subOrEmptyDict("integrationSchemes");
    
    TIntegrator_.reset
    (
        integrationScheme::New
        (
            "T",
            integrationSchemeDict_
        ).ptr()
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BedType>
inline Foam::ThermoBed<BedType>::ThermoBed
(
    const word& bedName,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo
)
:
    BedType(bedName, mesh, rho, U, thermo.thermo().mu(), g, thermo),
    thermo_(thermo),
    Tc_(thermo.thermo().T()),
    pc_(thermo.thermo().p()),
    Cpc_(thermo.thermo().Cp()),
    kappac_(thermo.thermo().kappa()),
    TMin_(this->constProperties_, 0.0),
    TMax_(this->constProperties_, vGreat),
    T0_(this->constProperties_, 0.0),
    Cp0_(this->constProperties_, 0.0),
    kp0_(this->constProperties_, 0.0),
    epsilon0_(this->constProperties_, 0.0),
    f0_(this->constProperties_, 0.0),
    bedT_
    (
        IOobject
        (
            IOobject::groupName(bedName, "T"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimTemperature, 0)
    ),
    bedCp_
    (
        IOobject
        (
            IOobject::groupName(bedName, "Cp"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Cp0())
    ),
    bedKp_
    (
        IOobject
        (
            IOobject::groupName(bedName, "kp"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimPower/dimLength/dimTemperature, kp0())
    ),
    heatTransferModel_(nullptr),
    TIntegrator_(nullptr),
    heatExplicit_(true),
    radiation_(false),
    radAreaP_(nullptr),
    radT4_(nullptr),
    radAreaPT4_(nullptr),
    hsTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                bedName + ":hsTrans",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy, 0)
        )
    ),
    hsCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                bedName + ":hsCoeff",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTemperature, 0)
        )
    )
{
    //- read scheme
    Istream& is = this->subModelProperties_.lookup("HeatSource");
    const word scheme(is);
    if (scheme == "semiImplicit")
    {
         heatExplicit_ = false;
    }
    else if (scheme == "explicit")
    {
        heatExplicit_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid scheme " << scheme << ". Valid schemes are "
            << "explicit and semiImplicit" << exit(FatalError);
    }
    
    setModels();
}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BedType>
void Foam::ThermoBed<BedType>::resetSourceTerms()
{
    hsTrans_->field() = 0.0;
    hsCoeff_->field() = 0.0;

    if (radiation_)
    {
        radAreaP_->field() = 0.0;
        radT4_->field() = 0.0;
        radAreaPT4_->field() = 0.0;
    }
}


template<class BedType>
void Foam::ThermoBed<BedType>::preSolve()
{
    
    Cpc() = thermo_.thermo().Cp();
    
    kappac() = thermo_.thermo().kappa();
    
    resetSourceTerms();
}


template<class BedType>
void Foam::ThermoBed<BedType>::calcSurfaceValues
(
    const label celli,
    const scalar T,
    scalar& Ts,
    scalar& rhos,
    scalar& mus,
    scalar& Pr,
    scalar& kappas
)
{
    // Surface temperature using two thirds rule
    Ts = (2.0*T + Tc()[celli])/3.0;

    if (Ts < TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting parcel surface temperature to "
                << TMin() <<  nl << endl;
        }

        Ts = TMin();
    }

    // Assuming thermo props vary linearly with T for small d(T)
    const scalar TRatio = Tc()[celli]/Ts;

    rhos = this->rho()[celli]*TRatio;

    mus = this->mu()[celli]/TRatio;
    kappas = kappac()[celli]/TRatio;

    Pr = Cpc()[celli]*mus/kappas;
    Pr = max(rootVSmall, Pr);
}
 
 
template<class BedType>
Foam::scalar Foam::ThermoBed<BedType>::calcHeatTransfer
(
    const label celli,
    const scalar dt,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW,
    const scalar Sh,
    const scalar Gc,
    scalar& dhsTrans,
    scalar& Sph
)
{
    const scalar d = this->dp()[celli];
    const scalar d2nd = this->dp2nd()[celli];
    const scalar As = this->areaS(d2nd);
    const scalar V = this->volume(d);
    const scalar rhop = this->rhop()[celli];
    const scalar m = V*rhop;
    const scalar T = bedT_[celli];
    const scalar cpbed = bedCp()[celli];

    // Calc heat transfer coefficient
    scalar htc = heatTransfer().htc(d, Re, Pr, kappa, NCpW);

    // Calculate the integration coefficients
    const scalar bcp = htc*As/(m*cpbed);
    const scalar acp = bcp*Tc()[celli];
    scalar ancp = Sh;
    if (radiation())
    {
        const scalar sigma = physicoChemical::sigma.value();

        ancp += As*epsilon0()*(Gc/4.0 - sigma*pow4(T));
    }
    ancp /= m*cpbed;

    // Integrate to find the new parcel temperature
    const scalar deltaT = TIntegrator().delta(T, dt, acp + ancp, bcp);
    const scalar deltaTncp = ancp*dt;
    const scalar deltaTcp = deltaT - deltaTncp;

    // Calculate the new temperature and the enthalpy transfer terms
    scalar Tnew = T + deltaT;
    Tnew = min(max(Tnew, TMin()), TMax());

    dhsTrans -= m*cpbed*deltaTcp;

    Sph = dt*m*cpbed*bcp;

    return Tnew;
}


template<class BedType>
void Foam::ThermoBed<BedType>::solveHeat(const scalar dt)
{
    preSolve();
    
    const labelList bedList = this->bedIDList();
    
    tmp<volScalarField> Gc
    (
        volScalarField::New
        (
            "G_coppy",
            this->mesh_,
            dimensionedScalar(dimPower/sqr(dimLength), 0)
        )
    ); 
    
    if (radiation())
    {
        Gc.ref() = this->mesh_.objectRegistry::template
                    lookupObject<volScalarField>("G");
    }
        

    forAll(bedList, i)
    {
        // Define local properties at beginning of time step
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        const label celli = bedList[i];
        const scalar np0 = this->particleNumber_[celli];

        // Store T for consistent radiation source
        const scalar T0 = bedT_[celli];


        // Calc surface values
        // ~~~~~~~~~~~~~~~~~~~
        scalar Ts, rhos, mus, Pr, kappas;
        calcSurfaceValues(celli, bedT_[celli], Ts, rhos, mus, Pr, kappas);

        // Reynolds number
        scalar Re = this->Re(rhos, this->U()[celli], this->dp2nd_[celli], mus);


        // Sources
        // ~~~~~~~

        // Explicit enthalpy source for particle
        scalar Sh = 0.0;

        // Linearised enthalpy source coefficient
        scalar Sph = 0.0;

        // Sensible enthalpy transfer from the particle to the carrier phase
        scalar dhsTrans = 0.0;


        // Heat transfer
        // ~~~~~~~~~~~~~
        
        // Radiation field 
        scalar G = 0.0;
        
        if (radiation())
        {
            G = Gc()[celli];
        }

        // Sum Ni*Cpi*Wi of emission species
        scalar NCpW = 0.0;
        
        

        // Calculate new particle temperature
        bedT_[celli] =
            this->calcHeatTransfer
            (
                celli,
                dt,
                Re,
                Pr,
                kappas,
                NCpW,
                Sh,
                G,
                dhsTrans,
                Sph
            );


        //  Accumulate carrier phase source terms
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Update sensible enthalpy transfer
        hsTrans()[celli] += np0*dhsTrans;

        // Update sensible enthalpy coefficient
        hsCoeff()[celli] += np0*Sph;

        // Update radiation fields
        if (radiation())
        {
            
            const scalar dp2nd = this->dp2nd_[celli];
            const scalar ap = this->areaP(dp2nd);
            const scalar T4 = pow4(T0);
            radAreaP()[celli] += dt*np0*ap;
            radT4()[celli] += dt*np0*T4;
            radAreaPT4()[celli] += dt*np0*ap*T4;
        }
        
    }
} 
 
 
 
 





// ************************************************************************* //
