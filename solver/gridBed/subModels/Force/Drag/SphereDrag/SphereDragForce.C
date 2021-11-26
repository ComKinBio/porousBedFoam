/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "SphereDragForce.H"

// * * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * //

template<class BedType>
Foam::scalar Foam::SphereDragForce<BedType>::CdRe(const scalar Re)
{
    if (Re > 1000.0)
    {
        return 0.424*Re;
    }
    else
    {
        return 24.0*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::SphereDragForce<BedType>::SphereDragForce
(
    BedType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<BedType>(owner, mesh, dict, typeName, false)
{}


template<class BedType>
Foam::SphereDragForce<BedType>::SphereDragForce
(
    const SphereDragForce<BedType>& df
)
:
    ParticleForce<BedType>(df)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::SphereDragForce<BedType>::~SphereDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
Foam::forceSuSp Foam::SphereDragForce<BedType>::calcCoupled
(
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
//     return forceSuSp(Zero, mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d())));
    return forceSuSp(Zero, mass*0.75*muc*CdRe(Re));
}


// ************************************************************************* //
