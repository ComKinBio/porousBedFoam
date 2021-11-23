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

#include "basicGravityBio2DBed.H"

namespace Foam
{
    defineTemplateTypeNameAndDebug(basicBio2DBed, 0);
}


namespace Foam
{
    defineTemplateTypeNameAndDebug(basicGravityBio2DBed, 0);
}

// Kinematic
#include "makeParcelMassTransferModels.H"
#include "makeBioParcelDryingModels.H"
#include "makeBioParcelDevolatilisationModels.H"
#include "makeBioParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Kinematic sub-models
makeParcelMassTransferModels(basicGravityBio2DBed);
makeBioParcelDryingModels(basicGravityBio2DBed);
makeBioParcelDevolatilisationModels(basicGravityBio2DBed);
makeBioParcelSurfaceReactionModels(basicGravityBio2DBed);
// ************************************************************************* //
