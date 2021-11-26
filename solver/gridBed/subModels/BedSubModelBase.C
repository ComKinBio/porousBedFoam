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

#include "BedSubModelBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::BedSubModelBase<BedType>::BedSubModelBase(BedType& owner)
:
    subModelBase(owner.outputProperties()),
    owner_(owner)
{}


template<class BedType>
Foam::BedSubModelBase<BedType>::BedSubModelBase
(
    BedType& owner,
    const dictionary& dict,
    const word& baseName,
    const word& modelType,
    const word& dictExt
)
:
    subModelBase
    (
        owner.outputProperties(),
        dict,
        baseName,
        modelType,
        dictExt
    ),
    owner_(owner)
{}


template<class BedType>
Foam::BedSubModelBase<BedType>::BedSubModelBase
(
    const word& modelName,
    BedType& owner,
    const dictionary& dict,
    const word& baseName,
    const word& modelType
)
:
    subModelBase
    (
        modelName,
        owner.outputProperties(),
        dict,
        baseName,
        modelType
    ),
    owner_(owner)
{}


template<class BedType>
Foam::BedSubModelBase<BedType>::BedSubModelBase
(
    const BedSubModelBase<BedType>& smb
)
:
    subModelBase(smb),
    owner_(smb.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::BedSubModelBase<BedType>::~BedSubModelBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
const BedType& Foam::BedSubModelBase<BedType>::owner() const
{
    return owner_;
}


template<class BedType>
BedType& Foam::BedSubModelBase<BedType>::owner()
{
    return owner_;
}


template<class BedType>
bool Foam::BedSubModelBase<BedType>::writeTime() const
{
    return
        active()
     && owner_.mesh().time().writeTime();
}


template<class BedType>
void Foam::BedSubModelBase<BedType>::write(Ostream& os) const
{
    writeEntry(os, "owner", owner_.bedName());
    subModelBase::write(os);
}


// ************************************************************************* //
