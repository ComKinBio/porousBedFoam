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

#include "DryingModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class BedType>
const Foam::wordList Foam::DryingModel<BedType>::
enthalpyTransferTypeNames
(
    IStringStream
    (
        "("
            "latentHeat "
            "enthalpyDifference"
        ")"
    )()
);


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BedType>
typename Foam::DryingModel<BedType>::enthalpyTransferType
Foam::DryingModel<BedType>::wordToEnthalpyTransfer(const word& etName)
const
{
    forAll(enthalpyTransferTypeNames, i)
    {
        if (etName == enthalpyTransferTypeNames[i])
        {
            return enthalpyTransferType(i);
        }
    }

    FatalErrorInFunction
        << "Unknown enthalpyType " << etName << ". Valid selections are:" << nl
        << enthalpyTransferTypeNames << exit(FatalError);

    return enthalpyTransferType(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::DryingModel<BedType>::DryingModel
(
    BedType& owner
)
:
    BedSubModelBase<BedType>(owner),
    enthalpyTransfer_(etLatentHeat),
    dMass_(0.0)
{}


template<class BedType>
Foam::DryingModel<BedType>::DryingModel
(
    const DryingModel<BedType>& pcm
)
:
    BedSubModelBase<BedType>(pcm),
    enthalpyTransfer_(pcm.enthalpyTransfer_),
    dMass_(pcm.dMass_)
{}


template<class BedType>
Foam::DryingModel<BedType>::DryingModel
(
    const dictionary& dict,
    BedType& owner,
    const word& type
)
:
    BedSubModelBase<BedType>(owner, dict, typeName, type),
    enthalpyTransfer_
    (
        wordToEnthalpyTransfer(this->coeffDict().lookup("enthalpyTransfer"))
    ),
    dMass_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::DryingModel<BedType>::~DryingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
const typename Foam::DryingModel<BedType>::enthalpyTransferType&
Foam::DryingModel<BedType>::enthalpyTransfer() const
{
    return enthalpyTransfer_;
}


template<class BedType>
Foam::scalar Foam::DryingModel<BedType>::dh
(
    const label idc,
    const label idl,
    const scalar p,
    const scalar T
) const
{
    return 0.0;
}


template<class BedType>
Foam::scalar Foam::DryingModel<BedType>::TMax
(
    const scalar p,
    const scalarField& X
) const
{
    return great;
}


template<class BedType>
Foam::scalar Foam::DryingModel<BedType>::Tvap(const scalarField& X) const
{
    return -great;
}


template<class BedType>
void Foam::DryingModel<BedType>::addToDryingMass(const scalar dMass)
{
    dMass_ += dMass;
}


template<class BedType>
void Foam::DryingModel<BedType>::info(Ostream& os)
{
    const scalar mass0 = this->template getBaseProperty<scalar>("mass");
    const scalar massTotal = mass0 + returnReduce(dMass_, sumOp<scalar>());

    Info<< "    Mass transfer phase change      = " << massTotal << nl;

    if (this->writeTime())
    {
        this->setBaseProperty("mass", massTotal);
        dMass_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DryingModelNew.C"

// ************************************************************************* //
