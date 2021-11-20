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

#include "bioBedFields.H"
#include "meshTools.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * *  Static Member  * * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(bioBedFields, 0);
}

template<>
const char* NamedEnum
<
    Foam::bioBedFields::particlePhase,
    4
>::names[] =
{
    "wet",
    "dry",
    "char",
    "ash"
};

const Foam::NamedEnum<Foam::bioBedFields::particlePhase, 4>
    Foam::bioBedFields::particlePhaseNames_;

    
const Foam::bioBedFields::particlePhase
Foam::bioBedFields::bedComponents[] = 
{
    Foam::bioBedFields::wet_, 
    Foam::bioBedFields::dry_, 
    Foam::bioBedFields::char_, 
    Foam::bioBedFields::ash_
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bioBedFields::bioBedFields
(
    const word& bedName,
    const fvMesh& mesh,
    const dictionary bioProperties_ 
)
:
    mesh_(mesh),
    bioProperties_(bioProperties_),
    moisture_(bioProperties_, "moisture"),
    ashContent_(bioProperties_, "ashContent"),
    ashFixedPorosity_(bioProperties_, "ashFixedPorosity"),
    rhopPtrList_(4),
    dp0PtrList_(4),
    T0PtrList_(4),
    Cp0PtrList_(4),
    kp0PtrList_(4),
    epsilon0PtrList_(4),
    f0PtrList_(4),
    particleNumberPtrList_(4),
    dpPtrList_(4),
    dp2ndPtrList_(4),
    massPtrList_(4),
    bedTPtrList_(4),
    bedCpPtrList_(4),
    bedKpPtrList_(4),
    conversionRatio_(4)
{
    for (const auto e : bedComponents)
    {
        word name(particlePhaseNames_[e]);
        
        word rhopname = name + "_rhop";
        rhopPtrList_[e].reset
        (
            new demandDrivenEntry<scalar>
            (
                bioProperties_,
                rhopname
            )
        );
        
        word dp0name = name + "_dp0";
        dp0PtrList_[e].reset
        (
            new demandDrivenEntry<scalar>
            (
                bioProperties_,
                dp0name
            )
        );
        
        word T0name = name + "_T0";
        T0PtrList_[e].reset
        (
            new demandDrivenEntry<scalar>
            (
                bioProperties_,
                T0name,
                0.0
            )
        );
        
        word Cp0name = name + "_Cp0";
        Cp0PtrList_[e].reset
        (
            new demandDrivenEntry<scalar>
            (
                bioProperties_,
                Cp0name,
                0.0
            )
        );
    
        word kp0name = name + "_kp0";
        kp0PtrList_[e].reset
        (
            new demandDrivenEntry<scalar>
            (
                bioProperties_,
                kp0name,
                0.0
            )
        );
    
        word epsilon0name = name + "_epsilon0";
        epsilon0PtrList_[e].reset
        (
            new demandDrivenEntry<scalar>
            (
                bioProperties_,
                epsilon0name,
                0.0
            )
        );
    
        word f0name = name + "_f0";
        f0PtrList_[e].reset
        (
            new demandDrivenEntry<scalar>
            (
                bioProperties_,
                f0name,
                0.0
            )
        );
    
        
        word particleNumbername = "particleNumber_" + name;
        particleNumberPtrList_[e].reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(bedName, particleNumbername),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, 0)
            )
        );
        
        word dpname = "dp_" + name;
        dpPtrList_[e].reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(bedName, dpname),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimLength, dp0e(e))
            )
        );
        
        word dp2ndname = "dp2nd_" + name;
        dp2ndPtrList_[e].reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(bedName, dp2ndname),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimLength, 0)
            )
        );
        
        word dpname = "mass_" + name;
        massPtrList_[e].reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(bedName, dpname),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (   dimMass, 
                    rhope(e)*particleNumbere(e)*pi/6.0*pow3(dpe(e))
                ) 
            )
        );
        
        
        word bedTname = "bedT_" + name;
        bedTPtrList_[e].reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(bedName, bedTname),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimTemperature, T0e(e))
            )
        );
        
        word bedCpname = "bedCp_" + name;
        bedCpPtrList_[e].reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(bedName, bedCpname),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    dimEnergy/dimMass/dimTemperature, 
                    Cp0e(e)
                )
            )
        );
        
        word bedKpname = "bedKp_" + name;
        bedKpPtrList_[e].reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(bedName, bedKpname),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (   dimPower/dimLength/dimTemperature, 
                    bedKpe(e)
                )
            )
        );
        
        word rationame = "registeredRatio_" + name;
        conversionRatio_[e].reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(bedName, rationame),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, 0)
            )
        );
    }
    
    ash_porosity() = dimensionedScalar(dimless, ashFixedPorosity_.value());
    dpe(ash_) = dpe(ash_)/cbrt(ash_porosity());
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




// ************************************************************************* //
