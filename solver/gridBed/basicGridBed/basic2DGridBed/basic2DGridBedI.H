/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fvmSup.H"
#include "integrationScheme.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //    
inline const Foam::word& 
Foam::basic2DGridBed::bedName() const
{
    return bedName_;
}

inline const Foam::fvMesh& 
Foam::basic2DGridBed::mesh() const
{
    return mesh_;
}


inline const Foam::IOdictionary& 
Foam::basic2DGridBed::bedProperties() const
{
    return bedProperties_;
}


inline const Foam::IOdictionary&
Foam::basic2DGridBed::outputProperties() const
{
    return outputProperties_;
}


inline Foam::IOdictionary& Foam::basic2DGridBed::outputProperties()
{
    return outputProperties_;
}


inline const Foam::dictionary& 
Foam::basic2DGridBed::particleProperties() const
{
    return particleProperties_;
}


inline Foam::Random& Foam::basic2DGridBed::rndGen() const
{
    return rndGen_;
}


inline const Foam::dictionary& 
Foam::basic2DGridBed::constProperties() const
{
    return constProperties_;
}


inline const Foam::dictionary& 
Foam::basic2DGridBed::subModelProperties() const
{
    return subModelProperties_;
}    


inline const Foam::volVectorField& 
Foam::basic2DGridBed::U() const
{
    return U_;
}


inline const Foam::volScalarField& 
Foam::basic2DGridBed::rho() const
{
    return rho_;
}


inline const Foam::volScalarField& 
Foam::basic2DGridBed::mu() const
{
    return mu_;
}


inline Foam::scalar
Foam::basic2DGridBed::alphaMin() const
{
    return alphaMin_;
}


inline Foam::scalar
Foam::basic2DGridBed::rhop0() const
{
    return rhop0_.value();
}


inline Foam::scalar
Foam::basic2DGridBed::sphericity() const
{
    return sphericity_.value();
}


inline const Foam::dimensionedVector& 
Foam::basic2DGridBed::g() const
{
    return g_;
}


// inline const Foam::basic2DGridBed::forceType&
inline const ParticleForce<Foam::basic2DGridBed>&
Foam::basic2DGridBed::forces() const
{
    return forces_();
}


inline Foam::ParticleForce<Foam::basic2DGridBed>&
Foam::basic2DGridBed::forces()
{
    return forces_();
}
  
    
inline bool Foam::basic2DGridBed::momentumExplicit() const
{
    return momentumExplicit_;
}


inline const Foam::integrationScheme&
Foam::basic2DGridBed::momentumIntegrator() const
{
    return momentumIntegrator_();
}


inline const Foam::tmp<Foam::volScalarField>
Foam::basic2DGridBed::bedVolume() const
{
    tmp<volScalarField> tbedVolume
    (
        volScalarField::New
        (
            bedName_ + ":bedVolume",
            mesh_,
            dimensionedScalar(dimVolume, 0)
        )
    );

    scalarField& bedVolume = tbedVolume.ref().primitiveFieldRef();
    forAll(bedIDList_, id)
    {
        label celli = bedIDList_[id];
        bedVolume[celli] = pi/6.0*pow3(dp_[celli])*particleNumber_[celli];
    }
    
    return tbedVolume;
}


inline const Foam::tmp<Foam::volScalarField>
Foam::basic2DGridBed::bedMass() const
{
    tmp<volScalarField> tbedMass
    (
        volScalarField::New
        (
            bedName_ + ":bedMass",
            mesh_,
            dimensionedScalar(dimMass, 0)
        )
    );

    scalarField& bedMass = tbedMass.ref().primitiveFieldRef();
    forAll(bedIDList_, id)
    {
        label celli = bedIDList_[id];
        bedMass[celli] = rhop_[celli]*pi/6.0*pow3(dp_[celli])*particleNumber_[celli];
    }
    
    return tbedMass;
}


inline const Foam::tmp<Foam::volScalarField>
Foam::basic2DGridBed::theta() const
{
    tmp<volScalarField> ttheta
    (
        volScalarField::New
        (
            bedName_ + ":theta",
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    scalarField& theta = ttheta.ref().primitiveFieldRef();
    forAll(bedIDList_, id)
    {
        label celli = bedIDList_[id];
        theta[celli] = rhop_[celli]*pi/6.0*pow3(dp_[celli])*particleNumber_[celli];
    }
    
    theta /= (mesh_.V()*rho_);
    
    return ttheta;
}


inline const Foam::tmp<Foam::volScalarField>
Foam::basic2DGridBed::beta() const
{
    tmp<volScalarField> tbeta
    (
        volScalarField::New
        (
            bedName_ + ":beta",
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    scalarField& beta = tbeta.ref().primitiveFieldRef();
    forAll(bedIDList_, id)
    {
        label celli = bedIDList_[id];
        beta[celli] =pi/6.0*pow3(dp_[celli])*particleNumber_[celli];
    }
    
    beta /= mesh_.V();
    
    return tbeta;
}

 
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::basic2DGridBed::UTrans()
{
    return UTrans_();
}


inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::basic2DGridBed::UTrans() const
{
    return UTrans_();
}


inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::basic2DGridBed::UCoeff()
{
    return UCoeff_();
}


inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::basic2DGridBed::UCoeff() const
{
    return UCoeff_();
}
        
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::basic2DGridBed::SU
(
    volVectorField& U, 
    const scalar dt
) const
{
    if (momentumExplicit_)
    {
        tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
        fvVectorMatrix& fvm = tfvm.ref();

        fvm.source() = -UTrans()/dt;

        return tfvm;
    }
    else
    {
        const volScalarField::Internal
            Vdt(mesh_.V()*dt);       

        return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
    }
}


inline Foam::tmp<Foam::fvVectorMatrix>
Foam::basic2DGridBed::SU
(
    volVectorField& U
) const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return SU(U, dt);
}


inline Foam::scalarList
Foam::basic2DGridBed::addedNum() const
{
    return addedNum_;
}


inline Foam::scalarList
Foam::basic2DGridBed::removedNum() const
{
    return absorbedNum_;
}                    


inline Foam::volScalarField
Foam::basic2DGridBed::particleNumber() const
{
    return particleNumber_;
}


inline Foam::volScalarField
Foam::basic2DGridBed::rhop() const
{
    return rhop_;
}


inline Foam::volScalarField
Foam::basic2DGridBed::dp() const
{
    return dp_;
}                    


inline Foam::volScalarField
Foam::basic2DGridBed::dp2nd() const
{
    return dp2nd_;
}


inline Foam::volScalarField
Foam::basic2DGridBed::alpha() const
{
    return alpha_;
} 


inline Foam::labelListList
Foam::basic2DGridBed::bedGridToFineGrid() const
{
    return bedGridToFineGrid_;
}


inline Foam::labelList 
Foam::basic2DGridBed::bedGridList() const
{
    return bedGridList_;
} 


inline Foam::scalarList
Foam::basic2DGridBed::bedGridVolume() const
{
    return bedGridVolume_;
}


inline Foam::labelListList
Foam::basic2DGridBed::bedGridColumn() const
{
    return bedGridColumn_;
}            


inline Foam::labelList
Foam::basic2DGridBed::bedGridBottom() const
{
    return bedGridBottom_;
}


inline Foam::labelList
Foam::basic2DGridBed::bedIDList() const
{
    return bedIDList_;
} 


inline Foam::labelList
Foam::basic2DGridBed::bedIDBottom() const
{
    return bedIDBottom_;
}


inline Foam::labelList 
Foam::basic2DGridBed::bedIDTop() const
{
    return bedIDTop_;
} 


inline Foam::labelList
Foam::basic2DGridBed::bedIDInternal() const
{
    return bedIDInternal_;
}


inline Foam::labelList
Foam::basic2DGridBed::bedIDTopBottom() const
{
    return bedIDTopBottom_;
}         
     

inline Foam::volScalarField& 
Foam::basic2DGridBed::particleNumber()
{
    return particleNumber_;
}


inline Foam::volScalarField& 
Foam::basic2DGridBed::rhop()
{
    return rhop_;
}


inline Foam::volScalarField& 
Foam::basic2DGridBed::dp()
{
    return dp_;
}


inline Foam::volScalarField& 
Foam::basic2DGridBed::dp2nd()
{
    return dp2nd_;
}


inline Foam::volScalarField& 
Foam::basic2DGridBed::alpha()
{
    return alpha_;
}


inline Foam::labelList& 
Foam::basic2DGridBed::bedIDList()
{
    return bedIDList_;
}                 


inline Foam::labelList& 
Foam::basic2DGridBed::bedIDBottom()
{
    return bedIDBottom_;
}               


inline Foam::labelList& 
Foam::basic2DGridBed::bedIDTop()
{
    return bedIDTop_;
}                        


inline Foam::labelList& 
Foam::basic2DGridBed::bedIDInternal()
{
    return bedIDInternal_;
}


inline Foam::labelList& 
Foam::basic2DGridBed::bedIDTopBottom()
{
    return bedIDTopBottom_;
}


inline void Foam::basic2DGridBed::alphaCalc()
{
    scalarField& alphaInter = alpha_.ref();
    scalarField dpInter = dp_.ref();
    scalarField Num = particleNumber_.ref();
    
    alphaInter = Num*pi/6.0*pow3(dpInter)/fluidGridCellVolumes_;
    
    alpha_ = max(1.0 - alpha_, alphaMin_);
}


inline Foam::scalar 
Foam::basic2DGridBed::mass(const label celli) const
{
    return rhop_[celli]*volume(celli);
}


inline Foam::scalar 
Foam::basic2DGridBed::volume(const label celli) const
{
    return volume(dp_[celli])*particleNumber_[celli];
}


inline Foam::scalar 
Foam::basic2DGridBed::volume(const scalar dp)
{
    return pi/6.0*pow3(dp);
}


inline Foam::scalar 
Foam::basic2DGridBed::areaP(const label celli) const
{
    return areaP(dp2nd_[celli])*particleNumber_[celli];
}


inline Foam::scalar 
Foam::basic2DGridBed::areaP(const scalar dp2nd)
{
    return 0.25*areaS(dp2nd);
}


inline Foam::scalar 
Foam::basic2DGridBed::areaS(const label celli) const
{
    return areaS(dp2nd_[celli])*particleNumber_[celli];
}


inline Foam::scalar 
Foam::basic2DGridBed::areaS(const scalar dp2nd)
{
    return pi*dp2nd*dp2nd;
}


inline Foam::scalar Foam::basic2DGridBed::Re
(
    const label celli
) const
{
    return Re(rho()[celli], U()[celli], dp2nd_[celli], mu()[celli]);
}


inline Foam::scalar Foam::basic2DGridBed::Re
(
    const scalar rhoc,
    const vector& Uc,
    const scalar d,
    const scalar muc
)
{
    return rhoc*mag(Uc)*d/max(muc, rootVSmall);
}
        
// ************************************************************************* //
