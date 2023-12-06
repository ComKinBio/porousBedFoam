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
#include "fvmLaplacian.H"
#include "fvmDdt.H"
#include "constants.H"

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
    
    this->subModelProperties().lookup("radiativeConduction") >> radiativeCond_;
    
    this->subModelProperties().lookup("heatConduction") >> heatConduction_;
    
    this->subModelProperties().lookup("cellBackground") >> cellBackground_;
    
    this->subModelProperties().lookup("Johnson") >> Johnson_;
    
    this->subModelProperties().lookup("ZBSModel") >> ZBS_;

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
    TMin_(this->constProperties_, "TMin", 0.0),
    TMax_(this->constProperties_, "TMax", vGreat),
    T0_(this->particleProperties_, "T0", 0.0),
    Cp0_(this->particleProperties_, "Cp0", 0.0),
    kp0_(this->particleProperties_, "kp0", 0.0),
    epsilon0_(this->particleProperties_, "epsilon0", 0.0),
    f0_(this->particleProperties_, "f0", 0.0),
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
    radiativeCond_(false),
    heatConduction_(false),
    cellBackground_(false),
    Johnson_(false),
    ZBS_(false),
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
    ),
    hsCond_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                bedName + ":hsCond",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime, 0)
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
    
    this->subModelProperties_.lookup("radiation") >> radiation_;

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
    if (radiation() && !radiativeCond())
    {
        const scalar sigma = physicoChemical::sigma.value();

        ancp += As*epsilon0()*(Gc/4.0 - sigma*pow4(T));
    }
    
    if (cellBackground())
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
void Foam::ThermoBed<BedType>::calcHeatConduction()
{
    const  scalar dt = this->mesh().time().deltaTValue();
   
    const volScalarField rhocp(this->rhop_*this->bedCp_);
    const volScalarField lambda(this->bedKp_);
    
    volScalarField Ts("TsNew", this->bedT_);
    volScalarField alphac("alphac", this->alpha_);    
    alphac = 1.0 - alphac;  
    
    volScalarField TT
    (
        IOobject
        (
            "TT",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    scalarField& TTInterFeildRef = TT.ref();
    scalarField& TsInterFeildRef = Ts.ref();
    
    TTInterFeildRef = TsInterFeildRef;
 
    if (ZBS())
    {
        //calculate bed properties
        volScalarField solidfraction
        (
            IOobject
            (
                "solidfraction",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar(dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        
        volScalarField ks
        (
            IOobject
            (
                "ks",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar(dimensionSet(1, 1, -3, -1, 0, 0, 0), 1e-8),
            zeroGradientFvPatchScalarField::typeName
        );

        //calculate fluid heat conductivity
        volScalarField kf
        (
            IOobject
            (
                "kf",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar(dimensionSet(1, 1, -3, -1, 0, 0, 0), 2.5143e-3),
            zeroGradientFvPatchScalarField::typeName
        );
        
        volScalarField dp
        (
            IOobject
            (
                "dpCond",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->dp_,
            zeroGradientFvPatchScalarField::typeName
        );
        
        forAll(solidfraction, i)
        {
            solidfraction[i] = alphac[i];
            if (solidfraction[i] < small)
            {
                solidfraction[i] = 1e-6;
            }
            
            if (dp[i] < small)
            {
                dp[i] = 0.001;
            }
            
            kf[i] = 2.5143e-3 + 7.7288e-5*Tc()[i] + 8.6248e-11*Foam::pow(Tc()[i], 2);
            ks[i] = max(lambda[i], 5.0286e-3); //two times kf
        }
        
        kf.correctBoundaryConditions();
        ks.correctBoundaryConditions();
        solidfraction.correctBoundaryConditions();
        dp.correctBoundaryConditions();

       //calculate bed effective conductivity  
        dimensionedScalar ems
        (
            "ems",
            dimensionSet(0, 0, 0, 0, 0, 0 ,0),
            0.85
        );
        
        dimensionedScalar Kphi
        (
            "Kphi",
            dimensionSet(0, 0, 0, 0, 0, 0 ,0),
            0.01
        );
        
        volScalarField Ts3
        (
            IOobject
            (
                "Ts3",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            Foam::pow(TT, 3),
            zeroGradientFvPatchScalarField::typeName
        );
        Ts3.correctBoundaryConditions();
        
        volScalarField k_eff_g = kf*(1. - Foam::sqrt(solidfraction));
       
        volScalarField kfs = kf/ks;

        volScalarField B = 1.25*Foam::pow((solidfraction/(1-solidfraction)), 10./9.);
        volScalarField B0 = 1.-kfs*B;
        volScalarField B1 = (1.-kfs)*B/Foam::sqr(B0)*Foam::log(1.0/kfs/B);
        volScalarField B2 = -(B+1.)/2;
        volScalarField B3 = -(B-1.)/B0;
//         volScalarField k_eff_c = kf*Foam::sqrt(solidfraction)*(2./B0)*(B1 + B2 + B3);
        
        volScalarField k_eff_c = kf*Foam::sqrt(solidfraction)*(2./B0)*(B1 + B2 + B3)*(1. - Kphi) + ks*Foam::sqrt(solidfraction)*Kphi;
               
        volScalarField k_eff_r_coeff = 4.*physicoChemical::sigma*this->dp_;

        volScalarField BBFE1 = (1. - Foam::sqrt(solidfraction))*(1. - solidfraction);
        volScalarField Lamda = ks/(4.*dp*physicoChemical::sigma*Ts3);
        volScalarField BBFE = BBFE1 + Foam::sqrt(solidfraction)/(2./ems - 1.)*(B + 1.)/B/(1. + 1/((2./ems - 1)*Lamda));
        volScalarField k_eff_r_coeffForT3 = BBFE*k_eff_r_coeff;
        
        if (Johnson())
        {
            scalar emm = ems.value();
            volScalarField FE =   5.08 - 24.63*solidfraction + 2.097*emm + 43.54*Foam::pow(solidfraction, 2) - 5.937*solidfraction*emm 
                        + 0.7391* Foam::pow(emm, 2) - 26.9*Foam::pow(solidfraction, 3) 
                        + 5.455*Foam::pow(solidfraction, 2)*emm - 0.6188*solidfraction*Foam::pow(emm, 2);
            
            k_eff_r_coeffForT3 = k_eff_r_coeff*FE;
        }
        
        volScalarField k_eff_r = k_eff_r_coeffForT3*Ts3;
        
        volScalarField ks_eff
        (
            IOobject
            (
                "ks_eff",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar(dimensionSet(1, 1, -3, -1, 0, 0, 0), 2.5143e-3),
            zeroGradientFvPatchScalarField::typeName
        );
        
        ks_eff = k_eff_g + k_eff_c + k_eff_r;
        ks_eff.correctBoundaryConditions();
        
        const scalar rhocpMax = max(rhocp).value();
        
        volScalarField rhocpM
        (
            IOobject
            (
                "rhocpM",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            rhocp,
            zeroGradientFvPatchScalarField::typeName
        );
        
        forAll(rhocpM,i)
        {
            if (rhocpM[i] < small)
            {
                rhocpM[i] = rhocpMax;
            }
        }
        
        rhocpM.correctBoundaryConditions();
        
        
        volScalarField DTs = ks_eff/rhocpM;
        
        solve
        (
            fvm::ddt(solidfraction, TT) 
        -
            fvm::laplacian
            (
                DTs, 
                TT
            )
        );

        
    }
    else
    {
        solve
        (
            fvm::ddt(alphac, rhocp, TT) 
        -
            fvm::laplacian
            (
                fvc::interpolate(alphac)
            *fvc::interpolate(lambda), 
                TT
            )
        );
    }
    
    TsInterFeildRef = TTInterFeildRef;

    tmp<volScalarField::Internal> thsCond
    (
        volScalarField::Internal::New
        (
            this->bedName() + ":thsCond",
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTime, 0)
        )
    );

    scalarField& thsCondRef = thsCond.ref();
    
    thsCondRef = rhocp*(Ts - this->bedT_)*this->mesh().V()*(1.0-this->alpha())/dt;   
   
    hsCond_() = thsCond;

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
