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

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::HeatTransferModel<BedType>::HeatTransferModel(BedType& owner)
:
    BedSubModelBase<BedType>(owner),
    BirdCorrection_(false)
{}


template<class BedType>
Foam::HeatTransferModel<BedType>::HeatTransferModel
(
    const dictionary& dict,
    BedType& owner,
    const word& type
)
:
    BedSubModelBase<BedType>(owner, dict, typeName, type),
    BirdCorrection_(this->coeffDict().lookup("BirdCorrection"))
{}


template<class BedType>
Foam::HeatTransferModel<BedType>::HeatTransferModel
(
    const HeatTransferModel<BedType>& htm
)
:
    BedSubModelBase<BedType>(htm),
    BirdCorrection_(htm.BirdCorrection_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::HeatTransferModel<BedType>::~HeatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
const Foam::Switch& Foam::HeatTransferModel<BedType>::BirdCorrection() const
{
    return BirdCorrection_;
}


template<class BedType>
Foam::scalar Foam::HeatTransferModel<BedType>::htc
(
    const scalar dp,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW
) const
{
    const scalar Nu = this->Nu(Re, Pr);

    scalar htc = Nu*kappa/dp;

    if (BirdCorrection_ && (mag(htc) > rootVSmall) && (mag(NCpW) > rootVSmall))
    {
        const scalar phit = min(NCpW/htc, 50);
        if (phit > 0.001)
        {
            htc *= phit/(exp(phit) - 1.0);
        }
    }

    return htc;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HeatTransferModelNew.C"

// ************************************************************************* //
