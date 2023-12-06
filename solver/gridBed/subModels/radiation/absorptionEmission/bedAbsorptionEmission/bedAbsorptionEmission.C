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

#include "bedAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(bed, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        bed,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::bed::bed
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    bedName_(coeffsDict_.lookup("bedName"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::bed::~bed()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::bed::aDisp(const label) const
{
    tmp<volScalarField> ta
    (
        volScalarField::New
        (
            "a",
            mesh_,
            dimensionedScalar(dimless/dimLength, 0)
        )
    );
    
    word fieldName = bedName_ + ".ap";
    
    ta.ref() += mesh_.objectRegistry::template
                            lookupObject<volScalarField>(fieldName);
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::bed::eDisp
(
    const label bandI
) const
{
    tmp<volScalarField> te
    (
        volScalarField::New
        (
            "e",
            mesh_,
            dimensionedScalar(dimless/dimLength, 0)
        )
    );

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::bed::EDisp
(
    const label bandI
) const
{
    tmp<volScalarField> tE
    (
        volScalarField::New
        (
            "E",
            mesh_,
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
        )
    );
    
    word fieldName = bedName_ + ".Ep";

    tE.ref() += mesh_.objectRegistry::template
                            lookupObject<volScalarField>(fieldName);
 
    // Total emission is 4 times the projected emission
    return 4*tE;
}


// ************************************************************************* //
