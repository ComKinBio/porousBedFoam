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

#include "BioBed.H"
#include "meshTools.H"

#include "bioBedFields.C"
#include "MassTransferModel.H"
#include "DryingModel.H"
#include "DevolatilisationModel.H"
#include "SurfaceReactionModel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //
template<class BedType>
void Foam::BioBed<BedType>::setModels()
{
    massTransferModel_.reset
    (
        MassTransferModel<BioBed<BedType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
    
    dryingModel_.reset
    (
        DryingModel<BioBed<BedType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
    
    devolatilisationModel_.reset
    (
        DevolatilisationModel<BioBed<BedType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
    
    surfaceReactionModel_.reset
    (
        SurfaceReactionModel<BioBed<BedType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class BedType>
void Foam::BioBed<BedType>::resetSourceTerms()
{
    forAll(rhoTrans_, i)
    {
        rhoTrans_[i].field() = 0.0;
    }
}


template<class BedType>
void Foam::BioBed<BedType>::preSolve()
{
    BedType::preSolve();
    
    resetSourceTerms();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::BioBed<BedType>::BioBed
(
    const word& bedName,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    const SLGThermo& thermo
)
:
    BedType(bedName, mesh, rho, U, g, thermo),
    bioBedFields(bedName, mesh, this->particleProperties_),
    TDevol_(this->constProperties_, "TDevol", 0.0),
    LDevol_(this->constProperties_, "LDevol", 0.0),
    wMinMass_(this->constProperties_, "massThreshold", 0.999),
    shrinkageFactorAlpha_(this->particleProperties_, "shrinkageFactorAlpha", 0.10),
    shrinkageFactorBeta_(this->particleProperties_, "shrinkageFactorBeta", 0.39),
    shrinkageFactorGamma_(this->particleProperties_, "shrinkageFactorGamma", 0.95),
    hRetentionCoeff_(this->constProperties_, "hRetentionCoeff", 0.0),
    massTransferModel_(nullptr),
    dryingModel_(nullptr),
    devolatilisationModel_(nullptr),
    surfaceReactionModel_(nullptr),
    massExplicit_(true),
    speciesExplicit_(true),
    randomConvesion_(false),
    gaussianConvesion_
    (
        readBool(this->subModelProperties_.lookup("gaussianDistributionConvesion"))
    ),
    rhoTrans_(thermo.carrier().species().size())
{
    Istream& isM = this->subModelProperties_.lookup("MassSource");
    const word schemeM(isM);
    if (schemeM == "semiImplicit")
    {
         massExplicit_ = false;
    }
    else if (schemeM == "explicit")
    {
        massExplicit_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid scheme " << schemeM << ". Valid schemes are "
            << "explicit and semiImplicit" << exit(FatalError);
    }

    Istream& isS = this->subModelProperties_.lookup("SpeciesSource");
    const word schemeS(isS);
    if (schemeS == "semiImplicit")
    {
         speciesExplicit_ = false;
    }
    else if (schemeS == "explicit")
    {
        speciesExplicit_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid scheme " << schemeS << ". Valid schemes are "
            << "explicit and semiImplicit" << exit(FatalError);
    }
    
    
    for (const auto e : bedComponents)
    {
        word name(particlePhaseNames_[e]);
        word dp2ndname = "dp2nd_" + name;

        wordList wrdList(2);
        wrdList[0] = mesh_.time().timeName();
        wrdList[1] = IOobject::groupName(this->bedName(), dp2ndname);

        fileName dp2ndnameFileName(wrdList);

        if (!exists(dp2ndnameFileName))
        {
            dp2nde(e) = dpe(e)/sqrt(this->sphericity());
        }
    }
  
  
    // Set storage for mass source fields and initialise to zero
    forAll(rhoTrans_, i)
    {
        const word& specieName = thermo.carrier().species()[i];
        rhoTrans_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    this->bedName() + ":rhoTrans_" + specieName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass, 0)
            )
        );
    }
    
    word moistureGaussianConv = "moistureGaussianConv.csv";
    
    wordList wordList(2);
    wordList[0] = mesh_.time().constant();
    wordList[1] = moistureGaussianConv;
    fileName moistureGaussianConvFileName(wordList);
    
    if (gaussianConvesion_)
    {
        moisturetable_ = interpolationTable<scalar>(moistureGaussianConvFileName);
        
        moisturetable_.check();
    } 
    
    setModels();
    
    updateBedThermo();
       
    componentUpdate();
    
    this->updateBedList();
    
    this->alphaCalc();
    
}

// * * * * * * * * * * * * * * * Main Function  * * * * * * * * * * * * * //

// solve thermo converison - void solveConversion
#include"BioBedMainFunc.H" 

template<class BedType>
void Foam::BioBed<BedType>::solveConversion()
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    if (this->radiativeCond() || this->heatConduction())
    {
        this->calcHeatConduction();
    }
    
    solveConversion(dt);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
void Foam::BioBed<BedType>::componentUpdate()
{
    updateBedThermo();

    tmp<volScalarField> tparticleNumber
    (
        volScalarField::New
        (
            this->bedName() + ":tparticleNumber",
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    scalarField& numberRef = tparticleNumber.ref();

    for (const auto e : bedComponents)
    {
        numberRef += particleNumbere(e);
    }
    
    //- update particle average diameter
    tmp<volScalarField> tdp
    (
        volScalarField::New
        (
            this->bedName() + ":tdp",
            this->mesh(),
            dimensionedScalar(dimLength, 0)
        )
    );

    scalarField& dpRef = tdp.ref();
    for (const auto e : bedComponents)
    {
        dpRef += pow3(dpe(e))*particleNumbere(e);
    }
    forAll(dpRef,i)
    {
        if (numberRef[i]>0)
        {
            dpRef[i] = cbrt(dpRef[i]/numberRef[i]);
        }
        else
        {
            dpRef[i] = small;
        }
    }
  
    //- update particle average diameter2
    tmp<volScalarField> tdp2
    (
        volScalarField::New
        (
            this->bedName() + ":tdp2",
            this->mesh(),
            dimensionedScalar(dimLength, 0)
        )
    );

    scalarField& dp2Ref = tdp2.ref();
    for (const auto e : bedComponents)
    {
        dp2Ref += sqr(dp2nde(e))*particleNumbere(e);
    }
    forAll(dp2Ref,i)
    {
        if (numberRef[i]>0)
        {
            dp2Ref[i] = sqrt(dp2Ref[i]/numberRef[i]);
        }
        else
        {
            dp2Ref[i] = small;
        }
    }
    
    //- update particle average rhop
    tmp<volScalarField> tmass
    (
        volScalarField::New
        (
            this->bedName() + ":tmass",
            this->mesh(),
            dimensionedScalar
            (
                dimMass, 
                0
            )
        )
    );

    scalarField& tmassRef = tmass.ref();
    for (const auto e : bedComponents)
    {
        tmassRef += masse(e);
    }
    
    tmp<volScalarField> trhop
    (
        volScalarField::New
        (
            this->bedName() + ":trhop",
            this->mesh(),
            dimensionedScalar
            (
                dimDensity, 
                0
            )
        )
    );
    
    scalarField& trhopRef = trhop.ref();
    forAll(trhopRef,i)
    {
        if (numberRef[i]>0)
        {
            trhopRef[i] = tmassRef[i]/numberRef[i]/this->volume(dpRef[i]);
        }
        else
        {
            trhopRef[i] = small;
        }
    }
    
    //- update particle average Cp
    tmp<volScalarField> tCp
    (
        volScalarField::New
        (
            this->bedName() + ":tCp",
            this->mesh(),
            dimensionedScalar
            (
                dimEnergy/dimMass/dimTemperature, 
                0
            )
        )
    );

    scalarField& tCpRef = tCp.ref();
    for (const auto e : bedComponents)
    {
        tCpRef += bedCpe(e)*masse(e);
    }
    forAll(tCpRef,i)
    {
        if (numberRef[i]>0)
        {
            tCpRef[i] = tCpRef[i]/tmassRef[i];
        }
        else
        {
            tCpRef[i] = small;
        }
    }
    
    //- update particle average Tp
    tmp<volScalarField> tTp
    (
        volScalarField::New
        (
            this->bedName() + ":tTp",
            this->mesh(),
            dimensionedScalar(dimTemperature, 0)
        )
    );

    scalarField& tTpRef = tTp.ref();
    for (const auto e : bedComponents)
    {
        tTpRef += bedTe(e)*masse(e)*bedCpe(e);
    }
    forAll(tTpRef,i)
    {
        if (numberRef[i]>0)
        {
            tTpRef[i] = tTpRef[i]/tmassRef[i]/tCpRef[i];
        }
        else
        {
            tTpRef[i] = this->Tc()[i];
        }
    }
    
    //- update particle average kp
    tmp<volScalarField> tkp
    (
        volScalarField::New
        (
            this->bedName() + ":tkp",
            this->mesh(),
            dimensionedScalar
            (
                dimPower/dimLength/dimTemperature, 
                0
            )
        )
    );
    
    const scalar sigma = physicoChemical::sigma.value();

    scalarField& tkpRef = tkp.ref();
    for (const auto e : bedComponents)
    {
        tkpRef += bedKpe(e)*masse(e);
    }
    forAll(tkpRef,i)
    {
        if (numberRef[i]>0)
        {
            tkpRef[i] = tkpRef[i]/tmassRef[i];
            
            if (this->radiativeCond_ && this->heatConduction_)
            {
                tkpRef[i] += 
                    4*this->alpha_[i]*sigma*this->epsilon0()*dpRef[i]*pow3(tTpRef[i]);
            }
            else if (this->radiativeCond_ && !this->heatConduction_)
            {
                tkpRef[i] = 
                    4*this->alpha_[i]*sigma*this->epsilon0()*dpRef[i]*pow3(tTpRef[i]);
            }
        }
        else
        {
            tkpRef[i] = small;
        }
    }
    
    
    this->particleNumber() = tparticleNumber;
    this->dp() = tdp;
    this->dp2nd() = tdp2;
    this->rhop() = trhop;
    this->bedCp() = tCp;
    this->bedKp() = tkp;
    this->bedT() = tTp;
    
    tmass.clear();
    
    this->alphaCalc();
}


template<class BedType>
void Foam::BioBed<BedType>::gaussianConvesionUpdate()
{
    // get ref fro the bed fields
    volScalarField& npWet = particleNumbere(wet_);
    volScalarField& npDry = particleNumbere(dry_);
    volScalarField& dWet = dpe(wet_);
    volScalarField& dDry = dpe(dry_);
    volScalarField& d2ndWet = dp2nde(wet_);
    volScalarField& d2ndDry = dp2nde(dry_);
    volScalarField& T0Wet = bedTe(wet_);
    volScalarField& T0Dry = bedTe(dry_);
    volScalarField& mass0Wet = masse(wet_);
    volScalarField& mass0Dry = masse(dry_);
    volScalarField& cpWet = bedCpe(wet_);
    volScalarField& cpDry = bedCpe(dry_);
    volScalarField& w = w_percent();
    volScalarField& gamma = gamma_percent();
    
    // loop check gaussianConvesion
    const labelList bedList = this->bedIDList();
    
    forAll(bedList, i)
    {
        const label celli = bedList[i];
        const scalar cellScale = cbrt(this->mesh().V()[celli]);
        const scalar numberThreshold = std::round(sqr(cellScale/dp0e(wet_)));
        
        if (npWet[celli] > numberThreshold)
        {
            const label celli = bedList[i];
            const scalar bedV = this->mesh().V()[celli]*this->alpha_[celli];
            scalar npi0 = bedV/this->volume(dp0e(wet_));
            npi0 = std::round(npi0);
            
            const scalar evaperated = (npi0*moisture() - npWet[celli]*w[celli])/npi0;
            
            if (evaperated > 0)
            {
                 const scalar dryRatio = moisturetable_(evaperated);
                 const scalar npWetGaussian = npi0 - std::round(npi0*dryRatio);

                 
                 if (0 < npWetGaussian < npWet[celli])
                 {
                    const scalar npWetToDry = npWet[celli] - npWetGaussian;
                    const scalar npWetNew = npWet[celli] - npWetToDry;
                                      
                    const scalar massWoodInWetNp1 = mass0Wet[celli]*(1-w[celli])/npWet[celli];
                    const scalar massWaterInWetNp1 = mass0Wet[celli]*w[celli]/npWet[celli];
                    const scalar transferedWood = massWoodInWetNp1*npWetToDry;
                    
                    const scalar newMoisture = mass0Wet[celli]*w[celli]/(mass0Wet[celli] - transferedWood);
                    const scalar newDryDp = dp0e(wet_)*cbrt(1-shrinkageFactorAlpha());
                    const scalar newDryDp2nd = newDryDp/sqrt(this->sphericity());
                     
                    T0Dry[celli] = (T0Dry[celli]*cpDry[celli]*mass0Dry[celli]
                                  + T0Wet[celli]*cpWet[celli]*transferedWood)\
                                  /(cpDry[celli]*mass0Dry[celli] + cpWet[celli]*transferedWood);
                    
                    dWet[celli] = dpAfterShrink(dWet[celli], mass0Wet[celli]/npWet[celli], \
                                  (massWaterInWetNp1*npWetToDry)/(moisture()*npWetNew), shrinkageFactorAlpha());
                    d2ndWet[celli] = dWet[celli]/sqrt(this->sphericity());

                    dDry[celli] = cbrt((pow3(dDry[celli])*npDry[celli] + pow3(newDryDp)*npWetToDry)\
                                 /(npDry[celli] + npWetToDry));
                    d2ndDry[celli] = sqrt((sqr(d2ndDry[celli])*npDry[celli] + sqr(newDryDp2nd)*npWetToDry)\
                                    /(npDry[celli] + npWetToDry));
                    
                    gamma[celli] = mass0Dry[celli]*gamma[celli]\
                                  /(mass0Dry[celli] + transferedWood);
                    mass0Dry[celli] += transferedWood;
                    
                    w[celli] = newMoisture;
                    mass0Wet[celli] -= transferedWood;
                    
                    npDry[celli] += npWetToDry;
                    npWet[celli] = npWetNew;
 
                 }
            }
        }
    }
    
    componentUpdate();
}


template<class BedType>
Foam::scalar Foam::BioBed<BedType>::calcHeatTransfer
(
    const scalar dt,
    const scalar d,
    const scalar d2nd,
    const scalar m,
    const scalar T,
    const scalar cpbed,
    const scalar Tc,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar shellR,
    const scalar NCpW,
    const scalar Sh,
    const scalar Gc,
    const scalar epsilon0,
    scalar& dhsTrans,
    scalar& Sph
)
{
    const scalar As = this->areaS(d2nd);

    // Calc heat transfer coefficient
    scalar htc = this->heatTransfer().htc(d, Re, Pr, kappa, NCpW);
    
    htc = htc/(shellR*As*htc + 1);

    // Calculate the integration coefficients
    const scalar bcp = htc*As/(m*cpbed);
    const scalar acp = bcp*Tc;
    scalar ancp = Sh;
    if (this->radiation() && !this->radiativeCond())
    {
        const scalar sigma = physicoChemical::sigma.value();

        ancp += As*epsilon0*(Gc/4.0 - sigma*pow4(T));
    }
    
    if (this->cellBackground())
    {
        const scalar sigma = physicoChemical::sigma.value();

        ancp += As*epsilon0*(Gc/4.0 - sigma*pow4(T));
    }
    
    ancp /= m*cpbed;

    // Integrate to find the new parcel temperature
    const scalar deltaT = this->TIntegrator().delta(T, dt, acp + ancp, bcp);
    const scalar deltaTncp = ancp*dt;
    const scalar deltaTcp = deltaT - deltaTncp;

    // Calculate the new temperature and the enthalpy transfer terms
    scalar Tnew = T + deltaT;
    Tnew = min(max(Tnew, this->TMin()), this->TMax());

    dhsTrans -= m*cpbed*deltaTcp;

    Sph = dt*m*cpbed*bcp;

    return Tnew;
}


template<class BedType>
void Foam::BioBed<BedType>::updateBedfieldsPrompt
(
    const label fromCell,
    const label toCell,
    const scalar number
)
{

    for (const auto e : bedComponents)
    {
        if (particleNumbere(e)[fromCell]>0)
        {
            scalar numbere = 
                number*particleNumbere(e)[fromCell]\
                /this->particleNumber()[fromCell];
            numbere = std::round(number);
                
            if (numbere > particleNumbere(e)[fromCell])
            {
                numbere = particleNumbere(e)[fromCell];
            }
            
            // Before transfer number records
            scalar npj = particleNumbere(e)[toCell];
            
            // Calc transfered mass
            scalar dmassij = numbere/particleNumbere(e)[fromCell]*masse(e)[fromCell];
            
            // update number
            particleNumbere(e)[fromCell] -= numbere;
            particleNumbere(e)[toCell] += numbere;    
            
            // update dp_, toCell only
            dpe(e)[toCell] = cbrt( (pow3(dpe(e)[toCell])*npj 
                                + pow3(dpe(e)[fromCell])*numbere)\
                                /particleNumbere(e)[toCell] );
            
            // update dp2nd_, toCell only
            dp2nde(e)[toCell] = sqrt( (sqr( dp2nde(e)[toCell])*npj 
                                    + sqr(dp2nde(e)[fromCell])*numbere)\
                                    /particleNumbere(e)[toCell] );
            
            // update conversion ratio
            if (e == wet_)
            {
                w_percent()[toCell] = ( dmassij*w_percent()[fromCell] 
                                    + masse(e)[toCell]*w_percent()[toCell] )\
                                    /(dmassij + masse(e)[toCell]);
            }
            else if (e == dry_)
            {
                gamma_percent()[toCell] = ( dmassij*gamma_percent()[fromCell] 
                                        + masse(e)[toCell]*gamma_percent()[toCell] )\
                                        /(dmassij + masse(e)[toCell]);
            }
            else if (e == char_)
            {
                eta_percent()[toCell] = ( dmassij*eta_percent()[fromCell] 
                                        + masse(e)[toCell]*eta_percent()[toCell] )\
                                        /(dmassij + masse(e)[toCell]);
            }
            
            // update temperature
            bedTe(e)[toCell] = ( bedTe(e)[fromCell]*dmassij 
                                + bedTe(e)[toCell]*masse(e)[toCell] )\
                                /(dmassij + masse(e)[toCell]);
            // update mass
            
            masse(e)[fromCell] -= dmassij;
            masse(e)[toCell] += dmassij;

            //- Constant, no update
            //- Bed specific heat capacity field [J/kg/K]
            //- and  heat conductivity field [J/s/m/K]        
            
            // update alpha in i,j
            this->alpha_[fromCell] = 
                max(1.0 - this->volume(fromCell)/this->mesh_.V()[fromCell], this->alphaMin_);
            this->alpha_[toCell] = 
                max(1.0 - this->volume(toCell)/this->mesh_.V()[toCell], this->alphaMin_);
        }
    }
    
    updateBedThermo();

}

template<class BedType>
void Foam::BioBed<BedType>::updateBedList()
{

    componentUpdate();
    
    BedType::updateBedList();
}
    
    

// template<class BedType>
// void Foam::BioBed<BedType>::updateBedfieldsPrompt
// (
//     const scalar number,
//     const label cell,
//     const scalar dp,
//     const scalar dp2nd
// )
// {
//     
// }
// 
// 
// template<class BedType>
// void Foam::BioBed<BedType>::updateBedfieldsPrompt
// (
//     const label cell,
//     const scalar number
// )
// {
//     
// }

// ************************************************************************* //
