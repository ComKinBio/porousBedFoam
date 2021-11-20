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

// * * * * * * * * * * * * * *  Static Member  * * * * * * * * * * * * * * //

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //
// for (const auto i : Components)
//     {
//         Info<< "particlePhase.i:"<<i<<nl<<endl;
//     }



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::BioBed<BedType>::BioBed
(
    const word& bedName,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo
)
:
    BedType(bedName, mesh, rho, U, g, thermo),
    bioBedFields(bedName, mesh, this->particleProperties_),
    TDevol_(this->constProperties_, 0.0),
    LDevol_(this->constProperties_, 0.0),
    shrinkageFactorAlpha_(this->particleProperties_, 0.10),
    shrinkageFactorBeta_(this->particleProperties_, 0.39),
    shrinkageFactorGamma_(this->particleProperties_, 0.95),
    hRetentionCoeff_(this->constProperties_, 0.0),
    massTransferModel_(nullptr),
    dryingModel_(nullptr),
    devolatilisationModel_(nullptr),
    surfaceReactionModel_(nullptr),
    massExplicit_(true),
    speciesExplicit_(true)
{
    for (const auto e : bedComponents)
    {
        dp2nde(e) = dpe(e)/sqrt(this->sphericity());
    }
    
    setModels();
}

// * * * * * * * * * * * * * * * Main Function  * * * * * * * * * * * * * //

// solve thermo converison - void solveConversion
#include"BioBedMainFunc.H" 

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
void Foam::BioBed<BedType>::dryingUpdate
(
    
)
{
    Random& rndGen = this->rndGen();
}


template<class BedType>
void Foam::BioBed<BedType>::devoUpdate()
{
    
}


template<class BedType>
void Foam::BioBed<BedType>::charBurnOutUpdate()
{
    
}


template<class BedType>
void Foam::BioBed<BedType>::componentUpdate()
{
    
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
    const scalar NCpW,
    const scalar Sh,
    const scalar Gc,
    const scalar epsilon0,
    scalar& dhsTrans,
    scalar& Sph
)
{
    const scalar As = this->areaS(d2nd);
    const scalar V = this->volume(d);

    // Calc heat transfer coefficient
    scalar htc = heatTransfer().htc(d, Re, Pr, kappa, NCpW);

    // Calculate the integration coefficients
    const scalar bcp = htc*As/(m*cpbed);
    const scalar acp = bcp*Tc;
    scalar ancp = Sh;
    if (this->radiation())
    {
        const scalar sigma = physicoChemical::sigma.value();

        ancp += As*epsilon0*(Gc/4.0 - sigma*pow4(T));
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


// ************************************************************************* //
