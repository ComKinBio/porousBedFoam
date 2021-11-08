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
inline Foam::MarkovBed<BedType>::MarkovBed
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
{}



// * * * * * * * MarkovBed Member Functions  * * * * * * * //

template<class BedType>
inline bool Foam::MarkovBed<BedType>::active() const
{
    return active_;
}


template<class BedType>
inline Foam::label Foam::MarkovBed<BedType>::typeId() const
{
    return typeId_;
}


template<class BedType>
inline Foam::scalar Foam::MarkovBed<BedType>::nParticle() const
{
    return nParticle_;
}


template<class BedType>
inline bool& Foam::MarkovBed<BedType>::active()
{
    return active_;
}


template<class BedType>
inline Foam::label& Foam::MarkovBed<BedType>::typeId()
{
    return typeId_;
}


template<class BedType>
inline Foam::scalar& Foam::MarkovBed<BedType>::nParticle()
{
    return nParticle_;
}


template<class BedType>
inline Foam::scalar Foam::MarkovBed<BedType>::volume() const
{
    return volume(d_);
}


template<class BedType>
inline Foam::scalar Foam::MarkovBed<BedType>::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}


// ************************************************************************* //
