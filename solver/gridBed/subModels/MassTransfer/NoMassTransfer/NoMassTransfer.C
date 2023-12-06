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

#include "NoMassTransfer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::NoMassTransfer<BedType>::NoMassTransfer
(
    const dictionary&,
    BedType& owner
)
:
    MassTransferModel<BedType>(owner)
{}


template<class BedType>
Foam::NoMassTransfer<BedType>::NoMassTransfer
(
    const NoMassTransfer<BedType>& htm
)
:
    MassTransferModel<BedType>(htm.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::NoMassTransfer<BedType>::~NoMassTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
bool Foam::NoMassTransfer<BedType>::active() const
{
    return false;
}


template<class BedType>
Foam::scalar Foam::NoMassTransfer<BedType>::Sh
(
    const scalar,
    const scalar
) const
{
    return 0.0;
}


template<class BedType>
Foam::scalar Foam::NoMassTransfer<BedType>::Sc() const
{
    return 1.0;
}


// ************************************************************************* //
