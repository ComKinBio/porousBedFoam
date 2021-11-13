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

// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

inline Foam::scalar
Foam::bioBedFields::rhope(particlePhase e) const
{
    return rhopPtrList_[e]().value();
}


inline Foam::scalar
Foam::bioBedFields::T0e(particlePhase e) const
{
    return T0PtrList_[e]().value();
}


inline Foam::scalar
Foam::bioBedFields::Cp0e(particlePhase e) const
{
    return Cp0PtrList_[e]().value();
}


inline Foam::scalar
Foam::bioBedFields::kp0e(particlePhase e) const
{
    return kp0PtrList_[e]().value();
}


inline Foam::scalar
Foam::bioBedFields::epsilon0e(particlePhase e) const
{
    return epsilon0PtrList_[e]().value();
}


inline Foam::scalar
Foam::bioBedFields::f0e(particlePhase e) const
{
    return f0PtrList_[e]().value();
}


inline Foam::volScalarField&
Foam::bioBedFields::particleNumbere(particlePhase e)
{
    return particleNumberPtrList_[e]();
}


inline Foam::volScalarField&
Foam::bioBedFields::dpe(particlePhase e)
{
    return dpPtrList_[e]();
}


inline Foam::volScalarField&
Foam::bioBedFields::dp2nde(particlePhase e)
{
    return dp2ndPtrList_[e]();
}


inline Foam::volScalarField&
Foam::bioBedFields::bedTe(particlePhase e)
{
    return bedTPtrList_[e]();
}


inline Foam::volScalarField&
Foam::bioBedFields::bedCpe(particlePhase e)
{
    return bedCpPtrList_[e]();
}


inline Foam::volScalarField&
Foam::bioBedFields::bedKpe(particlePhase e)
{
    return bedKpPtrList_[e]();
}


inline Foam::volScalarField&
Foam::bioBedFields::w_percent()
{
    return conversionRatio_[wet_]();
}


inline Foam::volScalarField&
Foam::bioBedFields::gamma_percent()
{
    return conversionRatio_[dry_]();
}


inline Foam::volScalarField&
Foam::bioBedFields::eta_percent()
{
    return conversionRatio_[char_]();
}

// ************************************************************************* //
