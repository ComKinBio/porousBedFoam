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

#include "NoDrying.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::NoDrying<BedType>::NoDrying
(
    const dictionary&,
    BedType& owner
)
:
    DryingModel<BedType>(owner)
{}


template<class BedType>
Foam::NoDrying<BedType>::NoDrying
(
    const NoDrying<BedType>& pcm
)
:
    DryingModel<BedType>(pcm.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::NoDrying<BedType>::~NoDrying()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
bool Foam::NoDrying<BedType>::active() const
{
    return false;
}


template<class BedType>
void Foam::NoDrying<BedType>::calculate
(
    const scalar dt,
    const label celli,
    const scalar Re,
    const scalar Nu,
    const scalar Sh,
    const scalar d,
    const scalar nu,
    const scalar T,
    const scalar Ts,
    const scalar pc,
    const scalar Tc,
    const scalarField& X,
    scalarField& dMassPC
) const
{
    // Nothing to do...
}


// ************************************************************************* //
