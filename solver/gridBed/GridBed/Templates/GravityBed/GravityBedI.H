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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BedType>
inline Foam::GravityBed<BedType>::GravityBed
(
    const word& bedName,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo
)
:
    BedType(bedName, mesh, rho, U, thermo.thermo().mu(), g, thermo),
    collapseSetting_(this->bedProperties_.subDict("collapseSetting")),
    alphaGCollapse_(collapseSetting_.lookupOrDefault<scalar>("collapseThreshold", 0.5)),
    alphaGCollapseMin_(collapseSetting_.lookupOrDefault<scalar>("collapseTo", 0.4)),
    alphaGCollapseTol_(collapseSetting_.lookupOrDefault<scalar>("collapseTolerance", 0.02)),
    skipCollapse_(false)
{}



// * * * * * * * GravityBed Member Functions  * * * * * * * //
template<class BedType>
inline const Foam::dictionary& 
Foam::GravityBed<BedType>::collapseSetting() const
{
    return collapseSetting_;
}


template<class BedType>
inline const Foam::scalar&
Foam::GravityBed<BedType>::alphaGCollapse() const
{
    return alphaGCollapse_;
}


template<class BedType>
inline const Foam::scalar&
Foam::GravityBed<BedType>::alphaGCollapseMin() const
{
    return alphaGCollapseMin_;
}


template<class BedType>
inline const Foam::scalar&
Foam::GravityBed<BedType>::alphaGCollapseTol() const
{
    return alphaGCollapseTol_;
}


// ************************************************************************* //
