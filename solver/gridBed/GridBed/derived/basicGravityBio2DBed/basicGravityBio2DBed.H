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

Class
    Foam::basicBio2DBed

Description
    Definition of basic biofuel bed with Gravity collapse in 2D

SourceFiles
    basicBio2DBed.H

\*---------------------------------------------------------------------------*/

#ifndef basicGravityBio2DBed_H
#define basicGravityBio2DBed_H

#include "contiguous.H"
#include "basic2DGridBed.H"
#include "ThermoBed.H"
#include "BioBed.H"
#include "GravityBed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
    typedef BioBed<ThermoBed<basic2DGridBed>> basicBio2DBed;

    template<>
    inline bool contiguous<basicBio2DBed>()
    {
        return true;
    }
    
    
    typedef GravityBed<BioBed<ThermoBed<basic2DGridBed>>> basicGravityBio2DBed;

    template<>
    inline bool contiguous<basicGravityBio2DBed>()
    {
        return true;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
