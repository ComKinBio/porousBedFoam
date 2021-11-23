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

Class
    Foam::bioBedFields

Description
    Class holds multi-components: wet wood, dry wood, char, and ash
    properties fields
    
SourceFiles
    bioBedFields.C


\*---------------------------------------------------------------------------*/

#ifndef bioBedFields_H
#define bioBedFields_H

#include "fvCFD.H"
#include "demandDrivenEntry.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bioBedFields Declaration
\*---------------------------------------------------------------------------*/


class bioBedFields
{

public:
        
    // Static enum declear
        
        enum particlePhase 
        {
            wet_, 
            dry_, 
            char_, 
            ash_
        };
         
        static const NamedEnum<particlePhase, 4> particlePhaseNames_;
            
        static const particlePhase bedComponents[];
    
    
protected:

    // Protected data

        //fluid Mesh class    
        const fvMesh& mesh_;

        //- Parcel properties dictionary
        const dictionary bioProperties_;
        
        //- Fields constant properties
        
            //- Initial moisture wet based
            demandDrivenEntry<scalar> moisture_;
            
            //- Initial ash content dry based
            demandDrivenEntry<scalar> ashContent_;
            
            //- Initial ash porosity is not changed
            demandDrivenEntry<scalar> ashFixedPorosity_;
        
        
        // bioField components properties, list size = 4
            
            //- Particle material density [kg/m^3]
            List<autoPtr<demandDrivenEntry<scalar>>> rhopPtrList_;
            
             //- Particle initial diameter
            List<autoPtr<demandDrivenEntry<scalar>>> dp0PtrList_;
            
            //- Particle initial temperature [K]
            List<autoPtr<demandDrivenEntry<scalar>>> T0PtrList_;

            //- Particle specific heat capacity [J/kg/K]
            List<autoPtr<demandDrivenEntry<scalar>>> Cp0PtrList_;
            
            //- Particle specific heat conductivity [J/s/m/K]
            List<autoPtr<demandDrivenEntry<scalar>>> kp0PtrList_;

            //- Particle emissivity [] (radiation)
            List<autoPtr<demandDrivenEntry<scalar>>> epsilon0PtrList_;

            //- Particle scattering factor [] (radiation)
            List<autoPtr<demandDrivenEntry<scalar>>> f0PtrList_;
            
       
        // bedComponents properties fields, list size = 4
            
            //- total components particle number
            List<autoPtr<volScalarField>> particleNumberPtrList_;

            //- volume average diameter of the particle fields [m]
            List<autoPtr<volScalarField>> dpPtrList_;
            
            //- surface average diameter of the particle fields [m]
            List<autoPtr<volScalarField>> dp2ndPtrList_;
            
            //- mass of the particle fields [m]
            List<autoPtr<volScalarField>> massPtrList_;
            
            //- Particle temperature field
            List<autoPtr<volScalarField>> bedTPtrList_;

            //- Bed specific heat capacity field [J/kg/K]
            List<autoPtr<volScalarField>> bedCpPtrList_;
            
            //- Bed specific heat conductivity field [J/s/m/K]
            List<autoPtr<volScalarField>> bedKpPtrList_;
            
            //- Particle registered properties
            List<autoPtr<volScalarField>> conversionRatio_;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("bioBedFields");


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        bioBedFields
        (
            const word& bedName,
            const fvMesh& mesh,
            const dictionary bioProperties_ 
        );
        
        //- Disallow default bitwise copy construction
        bioBedFields(const bioBedFields&) = delete;

        //- Destructor
        virtual ~bioBedFields()
        {}



    // Member Functions

        // Access
        
             //- Return const access to moisture wet based
            inline scalar moisture() const;
            
            //- Return const access to ash content dry based
            inline scalar ashContent() const;
            
            //- Return const access to initial ash porosity
            inline scalar ashFixedPorosity() const;
        
            //- Return const access to rhop [kg/m^3]
            inline scalar rhope(particlePhase e) const;
            
            //- Return const access to initial diameter [kg/m^3]
            inline scalar dp0e(particlePhase e) const;
            
            //- Return const access to initial temperature [K]
            inline scalar T0e(particlePhase e) const;
            
            //- Return const access to specific heat capacity [J/kg/K]
            inline scalar Cp0e(particlePhase e) const;
            
            //- Return const access to specific heat conductivity [J/s/m/K]
            inline scalar kp0e(particlePhase e) const;
            
            //- Return const access to particle emissivity [] (radiation)
            inline scalar epsilon0e(particlePhase e) const;
            
            //- Return const access to particle scattering factor [] (radiation)
            inline scalar f0e(particlePhase e) const;
            
            
        // Edit

            //- Return access to components particle number
            inline volScalarField& particleNumbere(particlePhase e);
            
            //- Return access to volume average diameter
            inline volScalarField& dpe(particlePhase e);
            
            //- Return access to particle surface average diameter
            inline volScalarField& dp2nde(particlePhase e);
            
            //- Return access to volume average diameter
            inline volScalarField& masse(particlePhase e);
            
            //- Return access to components temperature field
            inline volScalarField& bedTe(particlePhase e);
            
            //- Return access to components specific heat capacity field [J/kg/K]
            inline volScalarField& bedCpe(particlePhase e);
            
            //- Return access to components specific heat conductivity field [J/s/m/K]
            inline volScalarField& bedKpe(particlePhase e);
            
            //- Return access to water/dry_wood ratio in wet_
            inline volScalarField& w_percent();
            
            //- Return access to char/dry_wood ratio in dry_
            inline volScalarField& gamma_percent();
            
            //- Return access to ash/char ratio in char_
            inline volScalarField& eta_percent();
            
            //- Return access to ash porosity
            inline volScalarField& ash_porosity();

            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bioBedFieldsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
