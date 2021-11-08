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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
inline
Foam::BioBed<BedType>::constantProperties::constantProperties()
:
    dict_(dictionary::null),
    parcelTypeId_(dict_, -1),
    rhoMin_(dict_, 0.0),
    rho0_(dict_, 0.0),
    minParcelMass_(dict_, 0.0)
{}


template<class BedType>
inline Foam::BioBed<BedType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    dict_(cp.dict_),
    parcelTypeId_(cp.parcelTypeId_),
    rhoMin_(cp.rhoMin_),
    rho0_(cp.rho0_),
    minParcelMass_(cp.minParcelMass_)
{}


template<class BedType>
inline Foam::BioBed<BedType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    dict_(parentDict.subOrEmptyDict("constantProperties")),
    parcelTypeId_(dict_, "parcelTypeId", -1),
    rhoMin_(dict_, "rhoMin", 1e-15),
    rho0_(dict_, "rho0"),
    minParcelMass_(dict_, "minParcelMass", 1e-15)
{}


template<class BedType>
inline Foam::BioBed<BedType>::BioBed
(
    const polyMesh& mesh,
    const vector coarseGridOrigin_,
    const vector coarseGridVertex_,
    const scalar coarseGridSize_  
)
:
    BedType(mesh, coarseGridOrigin_, coarseGridVertex_, coarseGridSize_),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero)
{}




// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class BedType>
inline const Foam::dictionary&
Foam::BioBed<BedType>::constantProperties::dict() const
{
    return dict_;
}


template<class BedType>
inline Foam::label
Foam::BioBed<BedType>::constantProperties::parcelTypeId() const
{
    return parcelTypeId_.value();
}


template<class BedType>
inline Foam::scalar
Foam::BioBed<BedType>::constantProperties::rhoMin() const
{
    return rhoMin_.value();
}


template<class BedType>
inline Foam::scalar
Foam::BioBed<BedType>::constantProperties::rho0() const
{
    return rho0_.value();
}


template<class BedType>
inline Foam::scalar
Foam::BioBed<BedType>::constantProperties::minParcelMass() const
{
    return minParcelMass_.value();
}


// * * * * * * * BioBed Member Functions  * * * * * * * //

template<class BedType>
inline bool Foam::BioBed<BedType>::active() const
{
    return active_;
}


template<class BedType>
inline Foam::label Foam::BioBed<BedType>::typeId() const
{
    return typeId_;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::nParticle() const
{
    return nParticle_;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::d() const
{
    return d_;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::dTarget() const
{
    return dTarget_;
}


template<class BedType>
inline const Foam::vector& Foam::BioBed<BedType>::U() const
{
    return U_;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::rho() const
{
    return rho_;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::age() const
{
    return age_;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::tTurb() const
{
    return tTurb_;
}


template<class BedType>
inline const Foam::vector& Foam::BioBed<BedType>::UTurb() const
{
    return UTurb_;
}


template<class BedType>
inline bool& Foam::BioBed<BedType>::active()
{
    return active_;
}


template<class BedType>
inline Foam::label& Foam::BioBed<BedType>::typeId()
{
    return typeId_;
}


template<class BedType>
inline Foam::scalar& Foam::BioBed<BedType>::nParticle()
{
    return nParticle_;
}


template<class BedType>
inline Foam::scalar& Foam::BioBed<BedType>::d()
{
    return d_;
}


template<class BedType>
inline Foam::scalar& Foam::BioBed<BedType>::dTarget()
{
    return dTarget_;
}


template<class BedType>
inline Foam::vector& Foam::BioBed<BedType>::U()
{
    return U_;
}


template<class BedType>
inline Foam::scalar& Foam::BioBed<BedType>::rho()
{
    return rho_;
}


template<class BedType>
inline Foam::scalar& Foam::BioBed<BedType>::age()
{
    return age_;
}


template<class BedType>
inline Foam::scalar& Foam::BioBed<BedType>::tTurb()
{
    return tTurb_;
}


template<class BedType>
inline Foam::vector& Foam::BioBed<BedType>::UTurb()
{
    return UTurb_;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::mass() const
{
    return rho_*volume();
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::volume() const
{
    return volume(d_);
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::areaP() const
{
    return areaP(d_);
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::areaP(const scalar d)
{
    return 0.25*areaS(d);
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::areaS() const
{
    return areaS(d_);
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::areaS(const scalar d)
{
    return pi*d*d;
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::Re
(
    const scalar rhoc,
    const vector& U,
    const vector& Uc,
    const scalar d,
    const scalar muc
)
{
    return rhoc*mag(U - Uc)*d/max(muc, rootVSmall);
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::We
(
    const scalar rhoc,
    const vector& U,
    const vector& Uc,
    const scalar d,
    const scalar sigma
)
{
    return rhoc*magSqr(U - Uc)*d/max(sigma, rootVSmall);
}


template<class BedType>
inline Foam::scalar Foam::BioBed<BedType>::Eo
(
    const vector& g,
    const scalar rho,
    const scalar rhoc,
    const vector& U,
    const scalar d,
    const scalar sigma
)
{
    const vector dir = U/max(mag(U), rootVSmall);
    return mag(g & dir)*(rho - rhoc)*sqr(d)/max(sigma, rootVSmall);
}


// ************************************************************************* //
