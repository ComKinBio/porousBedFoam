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

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::ParticleForce<BedType>::ParticleForce
(
    BedType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType,
    const bool readCoeffs
)
:
    owner_(owner),
    mesh_(mesh),
    coeffs_
    (
        readCoeffs
      ? dict.optionalSubDict(forceType + "Coeffs")
      : dictionary::null
    )
{
    if (readCoeffs && coeffs_.isNull())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Force " << forceType << " must be specified as a dictionary"
            << exit(FatalIOError);
    }
}


template<class BedType>
Foam::ParticleForce<BedType>::ParticleForce(const ParticleForce& pf)
:
    owner_(pf.owner_),
    mesh_(pf.mesh_),
    coeffs_(pf.coeffs_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class BedType>
Foam::ParticleForce<BedType>::~ParticleForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
void Foam::ParticleForce<BedType>::cacheFields(const bool store)
{}


template<class BedType>
Foam::forceSuSp Foam::ParticleForce<BedType>::calcCoupled
(
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value;
    value.Su() = Zero;
    value.Sp() = 0.0;

    return value;
}


template<class BedType>
Foam::forceSuSp Foam::ParticleForce<BedType>::calcNonCoupled
(
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value;
    value.Su() = Zero;
    value.Sp() = 0.0;

    return value;
}


template<class BedType>
Foam::scalar Foam::ParticleForce<BedType>::massAdd
(
    const scalar mass
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleForceNew.C"

// ************************************************************************* //
