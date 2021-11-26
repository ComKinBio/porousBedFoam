/*---------------------------------------------------------------------------*\
 * =========                 |
 * \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 * \\    /   O peration     |
 *    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
 *    \\/     M anipulation  |
 * -------------------------------------------------------------------------------
 * License
 *    This file is part of OpenFOAM.
 * 
 *    OpenFOAM is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 * 
 *    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
 *    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *    for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * \*---------------------------------------------------------------------------*/

#include "StromNFDDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BedType>
Foam::StromNFDDevolatilisation<BedType>::
StromNFDDevolatilisation
(
    const dictionary& dict,
    BedType& owner
)
:
    DevolatilisationModel<BedType>(dict, owner, typeName),
    volgas_(this->coeffDict().lookup("volgas")),
    volatileToGasMap_(volgas_.size()),
    residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff"))),
    devolKinetic1_(this->coeffDict().lookup("devolKinetic1")),
    devolKinetic2_(this->coeffDict().lookup("devolKinetic2")),
    devolKinetic3_(this->coeffDict().lookup("devolKinetic3")),
    id_tar(owner.thermo().carrierId("tar"))
{
    if (volgas_.empty())
    {
        WarningIn
        (
            "Foam::StromNFDDevolatilisation<BedType>::"
            "StromNFDDevolatilisation"
            "("
            "const dictionary& dict, "
            "BedType& owner"
            ")"
        )   << "Devolatilisation model selected, but parameters are not well defined"
        << nl << endl;
    }
    else
    {
        
        Info << "Volatile info:" << endl;
        forAll(volgas_, i)
        {
            const word& specieName = volgas_[i].name();
            const label id = owner.thermo().carrierId(specieName);
            volatileToGasMap_[i] = id;
        }
       
        Info << " Kinetic data:" << endl;
        Info << "  wood -> gas A:" << devolKinetic1_[0] << " E:" << devolKinetic1_[1] <<endl;
        Info << "  wood -> tar A:" << devolKinetic2_[0] << " E:" << devolKinetic2_[1] <<endl;
        Info << "  wood -> char A:" << devolKinetic3_[0] << " E:" << devolKinetic3_[1] <<endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BedType>
Foam::StromNFDDevolatilisation<BedType>::
~StromNFDDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BedType>
void Foam::StromNFDDevolatilisation<BedType>::calculate
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
    bool done = true;
    
    // Model coefficients
    const scalar A1 = devolKinetic1_[0];
    const scalar E1 = devolKinetic1_[1];
    const scalar A2 = devolKinetic2_[0];
    const scalar E2 = devolKinetic2_[1];
    const scalar A3 = devolKinetic3_[0];
    const scalar E3 = devolKinetic3_[1];
    
    // Kinetic rate
    const scalar kappa1 = A1*exp(-E1/(RR*T));
    const scalar kappa2 = A2*exp(-E2/(RR*T));
    const scalar kappa3 = A3*exp(-E3/(RR*T));

    // initial calc
    const scalar massWood = mass;

    const scalar dmasswood_temp = dt*(kappa1+kappa2+kappa3)*massWood;   
    scalar tarprodmass = 0.0;
    
    if (dmasswood_temp < massWood /*&& massWood - dmasswood_temp >= massWood0*0.00005*/)
    {
        dMass = dmasswood_temp;
        tarprodmass = kappa2*massWood*dt;
        dMassChar = kappa3*massWood*dt;
        forAll(volgas_, i)
        {
            const label id_gas = volatileToGasMap_[i];
            dMassDV[id_gas] = volgas_[i].Y()*kappa1*massWood*dt;
        }
    }
    else
    {
        dMass = massWood;//-massWood0*0.00005;
        tarprodmass = kappa2*(dMass/(kappa1+kappa2+kappa3));    
        dMassChar = kappa3*(dMass/(kappa1+kappa2+kappa3));
        forAll(volgas_, i)
        {
            const label id_gas = volatileToGasMap_[i];
            dMassDV[id_gas] = volgas_[i].Y()*kappa1*(dMass/(kappa1+kappa2+kappa3));
        }    
    }

    dMassDV[id_tar] += tarprodmass;
   
    if (done && canCombust != -1)
    {
        canCombust = 1;
    }


}


// ************************************************************************* //
