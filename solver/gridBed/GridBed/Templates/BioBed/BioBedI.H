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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * BioBed Member Functions  * * * * * * * //

template<class BedType>
inline const Foam::SurfaceReactionModel<Foam::BioBed<BedType>>&
Foam::BioBed<BedType>::surfaceReaction() const
{
    return surfaceReactionModel_;
}


template<class BedType>
inline const Foam::MassTransferModel<Foam::BioBed<BedType>>&
Foam::BioBed<BedType>::massTransfer() const
{
    return massTransferModel_;
}


template<class BedType>
inline const Foam::DryingModel<Foam::BioBed<BedType>>&
Foam::BioBed<BedType>::drying() const
{
    return dryingModel_;
}


template<class BedType>
inline const Foam::DevolatilisationModel<Foam::BioBed<BedType>>&
Foam::BioBed<BedType>::devolatilisation() const
{
    return devolatilisationModel_;
}


template<class BedType>
inline Foam::SurfaceReactionModel<Foam::BioBed<BedType>>&
Foam::BioBed<BedType>::surfaceReaction()
{
    return surfaceReactionModel_();
}


template<class BedType>
inline Foam::MassTransferModel<Foam::BioBed<BedType>>&
Foam::BioBed<BedType>::massTransfer()
{
    return massTransferModel_();
}


template<class BedType>
inline Foam::DryingModel<Foam::BioBed<BedType>>&
Foam::BioBed<BedType>::drying()
{
    return dryingModel_();
}


template<class BedType>
inline Foam::DevolatilisationModel<Foam::BioBed<BedType>>&
Foam::BioBed<BedType>::devolatilisation()
{
    return devolatilisationModel_();
}


template<class BedType>
inline bool Foam::BioBed<BedType>::massExplicit() const
{
    return massExplicit_;
}


template<class BedType>
inline bool Foam::BioBed<BedType>::speciesExplicit() const
{
    return speciesExplicit_;
}


template<class BedType>
inline bool Foam::BioBed<BedType>::randomConvesion() const
{
    return randomConvesion_;
}


template<class BedType>
inline bool Foam::BioBed<BedType>::gaussianConvesion() const
{
    return gaussianConvesion_;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::TDevol() const
{
    return TDevol_.value();
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::LDevol() const
{
    return LDevol_.value();
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::wMinMass() const
{
    return wMinMass_.value();
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::shrinkageFactorAlpha() const
{
    return shrinkageFactorAlpha_.value();
}

template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::shrinkageFactorBeta() const
{
    return shrinkageFactorBeta_.value();
}

template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::shrinkageFactorGamma() const
{
    return shrinkageFactorGamma_.value();
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::hRetentionCoeff() const
{
    scalar value = hRetentionCoeff_.value();

    if ((value < 0) || (value > 1))
    {
        FatalErrorInFunction
            << "hRetentionCoeff must be in the range 0 to 1" << nl
            << exit(FatalError) << endl;
    }
    
    return value;
}


template<class BedType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::BioBed<BedType>::rhoTrans(const label i)
{
    return rhoTrans_[i];
}


template<class BedType>
inline
const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::BioBed<BedType>::rhoTrans() const
{
    return rhoTrans_;
}


template<class BedType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::BioBed<BedType>::rhoTrans()
{
    return rhoTrans_;
}


template<class BedType>
inline Foam::tmp<Foam::fvScalarMatrix> Foam::BioBed<BedType>::SYi
(
    const label i,
    volScalarField& Yi
) const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return SYi(i, Yi, dt);
}


template<class BedType>
inline Foam::tmp<Foam::fvScalarMatrix> Foam::BioBed<BedType>::SYi
(
    const label i,
    volScalarField& Yi,
    const scalar dt
) const
{
    if (speciesExplicit_)
    {
        tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Yi, dimMass/dimTime));
        fvScalarMatrix& fvm = tfvm.ref();

        fvm.source() = -rhoTrans_[i]/dt;

        return tfvm;
        
    }
    else
    {
        tmp<volScalarField> trhoTrans
        (
            volScalarField::New
            (
                this->bedName() + ":rhoTrans",
                this->mesh(),
                dimensionedScalar(dimMass/dimTime/dimVolume, 0)
            )
        );

        volScalarField& sourceField = trhoTrans.ref();
        
        dimensionedScalar dtdim(dimTime, dt);

        sourceField.primitiveFieldRef() =
            rhoTrans_[i]/(dtdim*this->mesh().V());

        const dimensionedScalar Yismall("Yismall", dimless, small);

        return
            fvm::Sp(neg(sourceField)*sourceField/(Yi + Yismall), Yi)
            + pos0(sourceField)*sourceField;
    }
}


template<class BedType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::BioBed<BedType>::Srho(const label i) const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return Srho(i, dt);
}


template<class BedType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::BioBed<BedType>::Srho(const label i, const scalar dt) const
{
    tmp<volScalarField::Internal> tRhoi
    (
        volScalarField::Internal::New
        (
            this->bedName() + ":rhoTrans",
            this->mesh(),
            dimensionedScalar
            (
                rhoTrans_[0].dimensions()/dimTime/dimVolume,
                0
            )
        )
    );

    scalarField& rhoi = tRhoi.ref();
    rhoi = rhoTrans_[i]/(dt*this->mesh().V());

    return tRhoi;
}


template<class BedType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::BioBed<BedType>::Srho() const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return Srho(dt);
}


template<class BedType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::BioBed<BedType>::Srho(const scalar dt) const
{
    tmp<volScalarField::Internal> trhoTrans
    (
        volScalarField::Internal::New
        (
            this->bedName() + ":rhoTrans",
            this->mesh(),
            dimensionedScalar
            (
                rhoTrans_[0].dimensions()/dimTime/dimVolume,
                0
            )
        )
    );

    scalarField& sourceField = trhoTrans.ref();
    forAll(rhoTrans_, i)
    {
        sourceField += rhoTrans_[i];
    }

    sourceField /= dt*this->mesh().V();
    
    return trhoTrans;

}


template<class BedType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::BioBed<BedType>::Srho(volScalarField& rho) const
{
    const scalar dt = this->mesh().time().deltaTValue();
    
    return Srho(rho, dt);
}


template<class BedType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::BioBed<BedType>::Srho(volScalarField& rho, const scalar dt) const
{
    tmp<volScalarField> trhoTrans
    (
        volScalarField::New
        (
            this->bedName() + ":rhoTrans",
            this->mesh(),
            dimensionedScalar(dimMass/dimTime/dimVolume, 0)
        )
    );

    scalarField& sourceField = trhoTrans.ref();

    if (massExplicit_)
    {

        tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
        fvScalarMatrix& fvm = tfvm.ref();

        forAll(rhoTrans_, i)
        {
            sourceField += rhoTrans_[i];
        }

        fvm.source() = -trhoTrans()/dt;

        return tfvm;
        
        
    }
    else
    {
        forAll(rhoTrans_, i)
        {
            sourceField += rhoTrans_[i];
        }
        
        dimensionedScalar dtdim(dimTime, dt);
        
        sourceField /= dtdim*this->mesh().V();

        return fvm::SuSp(trhoTrans()/rho, rho);
    }

}

template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::dpAfterShrink
(
    const scalar dp, 
    const scalar mass,
    const scalar dmass,
    const scalar shrinkage
) const
{
    const scalar VpOldcube = pow3(dp);

    const scalar rhopcube = mass/VpOldcube;
    
    const scalar dVpcube = dmass/rhopcube*shrinkage;
 
    return cbrt(VpOldcube - dVpcube);
}

// ************************************************************************* //
