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

#include "SingleKineticRateDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::SingleKineticRateDevolatilisation<BedType>::
SingleKineticRateDevolatilisation
(
    const dictionary& dict,
    BedType& owner
)
:
    DevolatilisationModel<BedType>(dict, owner, typeName),
    volatileData_(this->coeffDict().lookup("volatileData")),
    volatileToGasMap_(volatileData_.size()),
    residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff")))
{
    if (volatileData_.empty())
    {
        WarningInFunction
            << "Devolatilisation model selected, but no volatiles defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating volatile species:" << endl;

        // Determine mapping between active volatiles and bed gas components
        forAll(volatileData_, i)
        {
            const word& specieName = volatileData_[i].name();
            const label id = owner.thermo().carrierId(specieName);
            volatileToGasMap_[i] = id;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::SingleKineticRateDevolatilisation<BedType>::
~SingleKineticRateDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
void Foam::SingleKineticRateDevolatilisation<BedType>::calculate
(
    const scalar dt,
    const scalar age,
    const scalar mass0,
    const scalar mass,
    const scalar T,
    label& canCombust,
    scalar& dMass,
    scalarField& dMassDV,
    scalar& dMassChar
) const
{
    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];
        const scalar massVolatile = mass*volatileData_[i].Y();

        // Model coefficients
        const scalar A1 = volatileData_[i].A1();
        const scalar E = volatileData_[i].E();

        // Kinetic rate
        const scalar kappa = A1*exp(-E/(RR*T));

        // Mass transferred from particle to carrier gas phase
        dMassDV[id] = min(dt*kappa*massVolatile, massVolatile);
    }

    if (canCombust != -1)
    {
        canCombust = 1;
    }
}


// ************************************************************************* //
