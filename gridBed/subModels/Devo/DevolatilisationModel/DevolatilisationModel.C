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

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::DevolatilisationModel<BedType>::DevolatilisationModel
(
    BedType& owner
)
:
    BedSubModelBase<BedType>(owner),
    dMass_(0.0)
{}


template<class BedType>
Foam::DevolatilisationModel<BedType>::DevolatilisationModel
(
    const dictionary& dict,
    BedType& owner,
    const word& type
)
:
    BedSubModelBase<BedType>(owner, dict, typeName, type),
    dMass_(0.0)
{}


template<class BedType>
Foam::DevolatilisationModel<BedType>::DevolatilisationModel
(
    const DevolatilisationModel<BedType>& dm
)
:
    BedSubModelBase<BedType>(dm),
    dMass_(dm.dMass_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::DevolatilisationModel<BedType>::~DevolatilisationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
void Foam::DevolatilisationModel<BedType>::addToDevolatilisationMass
(
    const scalar dMass
)
{
    dMass_ += dMass;
}


template<class BedType>
void Foam::DevolatilisationModel<BedType>::info(Ostream& os)
{
    const scalar mass0 = this->template getBaseProperty<scalar>("mass");
    const scalar massTotal = mass0 + returnReduce(dMass_, sumOp<scalar>());

    Info<< "    Mass transfer devolatilisation  = " << massTotal << nl;

    if (this->writeTime())
    {
        this->setBaseProperty("mass", massTotal);
        dMass_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DevolatilisationModelNew.C"

// ************************************************************************* //
