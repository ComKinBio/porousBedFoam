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

#include "NoHeatTransfer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::NoHeatTransfer<BedType>::NoHeatTransfer
(
    const dictionary&,
    BedType& owner
)
:
    HeatTransferModel<BedType>(owner)
{}


template<class BedType>
Foam::NoHeatTransfer<BedType>::NoHeatTransfer
(
    const NoHeatTransfer<BedType>& htm
)
:
    HeatTransferModel<BedType>(htm.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::NoHeatTransfer<BedType>::~NoHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
bool Foam::NoHeatTransfer<BedType>::active() const
{
    return false;
}


template<class BedType>
Foam::scalar Foam::NoHeatTransfer<BedType>::Nu
(
    const scalar,
    const scalar
) const
{
    return 0.0;
}


template<class BedType>
Foam::scalar Foam::NoHeatTransfer<BedType>::Pr() const
{
    return 1.0;
}


// ************************************************************************* //
