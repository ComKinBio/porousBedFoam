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

#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * ThermoBed Member Functions  * * * * * * * //
template<class BedType>
inline const Foam::SLGThermo& Foam::ThermoBed<BedType>::thermo() const
{
    return thermo_;
}


template<class BedType>
inline const Foam::volScalarField& Foam::ThermoBed<BedType>::Tc() const
{
    return Tc_;
}


template<class BedType>
inline const Foam::volScalarField& Foam::ThermoBed<BedType>::pc() const
{
    return pc_;
}


template<class BedType>
inline Foam::volScalarField Foam::ThermoBed<BedType>::Cpc() const
{
    return Cpc_;
}


template<class BedType>
inline Foam::volScalarField Foam::ThermoBed<BedType>::kappac() const
{
    return kappac_;
}


template<class BedType>
inline Foam::scalar
Foam::ThermoBed<BedType>::T0() const
{
    return T0_.value();
}


template<class BedType>
inline Foam::scalar
Foam::ThermoBed<BedType>::TMin() const
{
    return TMin_.value();
}


template<class BedType>
inline Foam::scalar
Foam::ThermoBed<BedType>::TMax() const
{
    return TMax_.value();
}


template<class BedType>
inline Foam::scalar
Foam::ThermoBed<BedType>::Cp0() const
{
    return Cp0_.value();
}


template<class BedType>
inline Foam::scalar
Foam::ThermoBed<BedType>::kp0() const
{
    return kp0_.value();
}


template<class BedType>
inline Foam::scalar
Foam::ThermoBed<BedType>::epsilon0() const
{
    return epsilon0_.value();
}


template<class BedType>
inline Foam::scalar
Foam::ThermoBed<BedType>::f0() const
{
    return f0_.value();
}


template<class BedType>
inline Foam::volScalarField Foam::ThermoBed<BedType>::bedT() const
{
    return bedT_;
}


template<class BedType>
inline Foam::volScalarField Foam::ThermoBed<BedType>::bedCp() const
{
    return bedCp_;
}


template<class BedType>
inline Foam::volScalarField Foam::ThermoBed<BedType>::bedKp() const
{
    return bedKp_;
}

template<class BedType>
inline const Foam::HeatTransferModel<Foam::ThermoBed<BedType>>&
Foam::ThermoBed<BedType>::heatTransfer() const
{
    return heatTransferModel_;
}

template<class BedType>
inline bool Foam::ThermoBed<BedType>::heatExplicit() const
{
    return heatExplicit_;
}

template<class BedType>
inline const Foam::integrationScheme&
Foam::ThermoBed<BedType>::TIntegrator() const
{
    return TIntegrator_();
}


template<class BedType>
inline bool Foam::ThermoBed<BedType>::radiation() const
{
    return radiation_;
}


template<class BedType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::radAreaP()
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaP_();
}


template<class BedType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::radAreaP() const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaP_();
}


template<class BedType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::radT4()
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radT4_();
}


template<class BedType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::radT4() const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radT4_();
}


template<class BedType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::radAreaPT4()
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPT4_();
}


template<class BedType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::radAreaPT4() const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return radAreaPT4_();
}


template<class BedType>
inline Foam::volScalarField& Foam::ThermoBed<BedType>::Cpc()
{
    return Cpc_;
}        


template<class BedType>
inline Foam::volScalarField& Foam::ThermoBed<BedType>::kappac()
{
    return kappac_;
} 


template<class BedType>
inline Foam::volScalarField& Foam::ThermoBed<BedType>::bedT()
{
    return bedT_;
} 


template<class BedType>
inline Foam::volScalarField& Foam::ThermoBed<BedType>::bedCp()
{
    return bedCp_;
} 


template<class BedType>
inline Foam::volScalarField& Foam::ThermoBed<BedType>::bedKp()
{
    return bedKp_;
} 


template<class BedType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::hsTrans()
{
    return hsTrans_();
}


template<class BedType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::hsTrans() const
{
    return hsTrans_();
}


template<class BedType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::hsCoeff()
{
    return hsCoeff_();
}


template<class BedType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoBed<BedType>::hsCoeff() const
{
    return hsCoeff_();
}


template<class BedType>
inline Foam::tmp<Foam::fvScalarMatrix> 
Foam::ThermoBed<BedType>::Sh(volScalarField& hs) const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return Sh(hs, dt);
}


template<class BedType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ThermoBed<BedType>::Sh(volScalarField& hs, const scalar dt) const
{
    if (heatExplicit_)
    {
        tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(hs, dimEnergy/dimTime));
        fvScalarMatrix& fvm = tfvm.ref();

        fvm.source() = -hsTrans()/dt;

        return tfvm;
    }
    else
    {        
        dimensionedScalar dtdim(dimTime, dt);
        
        const volScalarField Cp(thermo_.thermo().Cp());
        const volScalarField::Internal
            Vdt(this->mesh().V()*dtdim);

        return
            hsTrans()/Vdt
            - fvm::SuSp(hsCoeff()/(Cp*Vdt), hs)
            + hsCoeff()/(Cp*Vdt)*hs;
    }
}


template<class BedType>
inline Foam::tmp<Foam::volScalarField> 
Foam::ThermoBed<BedType>::Ep() const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return Ep(dt);
}


template<class BedType>
inline Foam::tmp<Foam::volScalarField> 
Foam::ThermoBed<BedType>::Ep(const scalar dt) const
{
    tmp<volScalarField> tEp
    (
        volScalarField::New
        (
            this->bedName() + ":radiation:Ep",
            this->mesh(),
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
        )
    );

    if (radiation_)
    {
        scalarField& Ep = tEp.ref().primitiveFieldRef();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = epsilon0();
        const scalarField& sumAreaPT4 = radAreaPT4_->field();

        Ep = sumAreaPT4*epsilon*physicoChemical::sigma.value()/V/dt;
    }

    return tEp;
}


template<class BedType>
inline Foam::tmp<Foam::volScalarField> 
Foam::ThermoBed<BedType>::ap() const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return ap(dt);
}


template<class BedType>
inline Foam::tmp<Foam::volScalarField> 
Foam::ThermoBed<BedType>::ap(const scalar dt) const
{
    tmp<volScalarField> tap
    (
        volScalarField::New
        (
            this->bedName() + ":radiation:ap",
            this->mesh(),
            dimensionedScalar(dimless/dimLength, 0)
        )
    );

    if (radiation_)
    {
        scalarField& ap = tap.ref().primitiveFieldRef();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = epsilon0();
        const scalarField& sumAreaP = radAreaP_->field();

        ap = sumAreaP*epsilon/V/dt;
    }

    return tap;
}


template<class BedType>
inline Foam::tmp<Foam::volScalarField> 
Foam::ThermoBed<BedType>::sigmap() const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return sigmap(dt);
}


template<class BedType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoBed<BedType>::sigmap(const scalar dt) const
{
    tmp<volScalarField> tsigmap
    (
        volScalarField::New
        (
            this->bedName() + ":radiation:sigmap",
            this->mesh(),
            dimensionedScalar(dimless/dimLength, 0)
        )
    );

    if (radiation_)
    {
        scalarField& sigmap = tsigmap.ref().primitiveFieldRef();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = epsilon0();
        const scalar f = f0();
        const scalarField& sumAreaP = radAreaP_->field();

        sigmap = sumAreaP*(1.0 - f)*(1.0 - epsilon)/V/dt; 
    }

    return tsigmap;
}




// ************************************************************************* //
