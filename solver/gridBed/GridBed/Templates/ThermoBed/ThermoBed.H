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
    Foam::ThermoBed

Description
    ThermoBed class with with thermol properties and heat transfer models

    Sub-models include:
    - heat transfer

SourceFiles
    ThermoBedI.H
    ThermoBed.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoBed_H
#define ThermoBed_H

#include "IOstream.H"
#include "autoPtr.H"
#include "demandDrivenEntry.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class BedType>
class ThermoBed;


template<class bedType>
class HeatTransferModel;

/*---------------------------------------------------------------------------*\
                         Class ThermoBed Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class ThermoBed
:
    public BedType
{

public:

    // Public Typedefs

        //- Type of bed this bed was instantiated for
        typedef BedType bedType;

        //- Convenience typedef for this bed type
        typedef ThermoBed<BedType> thermoBedType;
        

protected:

    // Protected data
   
        // References to the carrier gas fields

            //- SLG thermodynamics package
            const SLGThermo& thermo_;

            //- Temperature [K]
            const volScalarField& Tc_;

            //- Pressure [Pa]
            const volScalarField& pc_;
            
            //- Update in presolve
            //- Local copy of carrier specific heat field
            //  Cp not stored on carrier thermo, but returned as tmp<...>
            volScalarField Cpc_;

            //- Local copy of carrier thermal conductivity field
            //  kappa not stored on carrier thermo, but returned as tmp<...>
            volScalarField kappac_;
        

        // bed constant properties
            
            //- Minimum temperature [K]
            demandDrivenEntry<scalar> TMin_;

            //- Maximum temperature [K]
            demandDrivenEntry<scalar> TMax_;
       
            //- Particle initial temperature [K]
            demandDrivenEntry<scalar> T0_;

            //- Particle specific heat capacity [J/kg/K]
            demandDrivenEntry<scalar> Cp0_;
            
            //- Particle specific heat conductivity [J/s/m/K]
            demandDrivenEntry<scalar> kp0_;

            //- Particle emissivity [] (radiation)
            demandDrivenEntry<scalar> epsilon0_;

            //- Particle scattering factor [] (radiation)
            demandDrivenEntry<scalar> f0_;
            
            
        // bed properties fields
                    
            //- Particle temperature field
            volScalarField bedT_;

            //- Bed specific heat capacity field [J/kg/K]
            volScalarField bedCp_;
            
            //- Bed specific heat conductivity field [J/s/m/K]
            volScalarField bedKp_;
         
        
        // bed submodels

            //- Heat transfer model
            autoPtr<HeatTransferModel<ThermoBed<BedType>>>
                heatTransferModel_;
            
         // Reference to the particle integration schemes

            //- Temperature integration
            autoPtr<integrationScheme> TIntegrator_;
            
            //- Heat source term shceme
            bool heatExplicit_;
            
            
        // Modelling options

            //- Include radiation
            Switch radiation_;

            //- Radiation sum of parcel projected areas
            autoPtr<volScalarField::Internal> radAreaP_;

            //- Radiation sum of parcel temperature^4
            autoPtr<volScalarField::Internal> radT4_;

            //- Radiation sum of parcel projected areas * temperature^4
            autoPtr<volScalarField::Internal> radAreaPT4_;


        // Sources

            //- Sensible enthalpy transfer [J/kg]
            autoPtr<volScalarField::Internal> hsTrans_;

            //- Coefficient for carrier phase hs equation [W/K]
            autoPtr<volScalarField::Internal> hsCoeff_;


    // Protected Member Functions

         // Initialisation

            //- Set cloud sub-models
            void setModels();
            
        
        // bed calculate functions

            //- PreSolve 
            void preSolve();
            
            //- resetSourceTerms 
            void resetSourceTerms();

        
       

public:

    // Static Data Members

        //- Runtime type information
        TypeName("ThermoBed");


    // Constructors

        //- Construct from mesh, coordinates and topology        
        //  Other properties initialised as null
        ThermoBed
        (
            const word& bedName,
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo
        );

//         //- Construct as a copy
//         ThermoBed(const ThermoBed& p);
// 
//         //- Construct as a copy
//         ThermoBed(const ThermoBed& p, const polyMesh& mesh);

        //- Destructor
        virtual ~ThermoBed()
        {}



    // Member Functions

        // Access
            
            //- Return const access to thermo package
            inline const SLGThermo& thermo() const;

            //- Return const access to the carrier temperature field
            inline const volScalarField& Tc() const;

            //- Return const access to the carrier prressure field
            inline const volScalarField& pc() const;
            
            //- Return const access to carrier specific heat field
            inline volScalarField Cpc() const;
            
            //- Return const access to carrier thermal conductivity field
            inline volScalarField kappac() const;
            
            //- Return const access to minimum temperature [K]
            inline scalar TMin() const;

            //- Return const access to maximum temperature [K]
            inline scalar TMax() const;
            
             //- Return const access to the particle initial temperature [K]
            inline scalar T0() const;

            //- Return const access to the particle specific heat capacity
            //  [J/kg/K]
            inline scalar Cp0() const;
            
            //- Return const access to the particle specific heat conductivity
            //  [J/s/m/K]
            inline scalar kp0() const;

            //- Return const access to the particle emissivity []
            //  Active for radiation only
            inline scalar epsilon0() const;

            //- Return const access to the particle scattering factor []
            //  Active for radiation only
            inline scalar f0() const;
            
            //- Return const access to bed temperature field
            inline volScalarField bedT() const;
            
            //- Return const access to bed specific heat capacity field
            inline volScalarField bedCp() const;
            
            //- Return const access to bed specific heat conductivity field
            inline volScalarField bedKp() const;
            
        
        // Sub-models
            
            //-Return reference to temperature integration
            inline const integrationScheme& TIntegrator() const; 
           
            //- Heat source scheme
            inline bool heatExplicit() const;

            //- Return reference to heat transfer model
            inline const HeatTransferModel<ThermoBed<bedType>>&
                    heatTransfer() const;
                    
            //- Radiation flag
            inline bool radiation() const;

            //- Radiation sum of parcel projected areas [m^2]
            inline const volScalarField::Internal&
                radAreaP() const;

            //- Radiation sum of parcel temperature^4 [K4]
            inline const volScalarField::Internal& radT4() const;

            //- Radiation sum of parcel temperature^4 [m2K4]
            inline const volScalarField::Internal&
                radAreaPT4() const;
                
            
        // Edit
                
            //- Return const access to carrier specific heat field
            inline volScalarField& Cpc();
            
            //- Return const access to carrier thermal conductivity field
            inline volScalarField& kappac();
            
            //- Return const access to bed temperature field
            inline volScalarField& bedT();
            
            //- Return const access to bed specific heat capacity field
            inline volScalarField& bedCp();
            
            //- Return const access to bed specific heat conductivity field
            inline volScalarField& bedKp();

            //- Radiation sum of parcel projected areas [m^2]
            inline volScalarField::Internal& radAreaP();

            //- Radiation sum of parcel temperature^4 [K4]
            inline volScalarField::Internal& radT4();
            
            //- Radiation sum of parcel projected area*temperature^4 [m2K4]
            inline volScalarField::Internal& radAreaPT4();


        // Sources

            // Enthalpy

                //- Sensible enthalpy transfer [J/kg]
                inline volScalarField::Internal& hsTrans();

                //- Sensible enthalpy transfer [J/kg]
                inline const volScalarField::Internal&
                    hsTrans() const;

                //- Return coefficient for carrier phase hs equation
                inline volScalarField::Internal& hsCoeff();

                //- Return const coefficient for carrier phase hs equation
                inline const volScalarField::Internal&
                    hsCoeff() const;

                //- Return sensible enthalpy source term [J/kg/m^3/s]
                inline tmp<fvScalarMatrix> Sh(volScalarField& hs, const scalar dt) const;
                
                //- Return sensible enthalpy source term [J/kg/m^3/s]
                inline tmp<fvScalarMatrix> Sh(volScalarField& hs) const;


            // Radiation

                //- Return tmp equivalent particulate emission
                inline tmp<volScalarField> Ep() const;

                //- Return tmp equivalent particulate absorption
                inline tmp<volScalarField> ap() const;

                //- Return tmp equivalent particulate scattering factor
                inline tmp<volScalarField> sigmap() const;
                
                //- Return tmp equivalent particulate emission
                inline tmp<volScalarField> Ep(const scalar dt) const;

                //- Return tmp equivalent particulate absorption
                inline tmp<volScalarField> ap(const scalar dt) const;

                //- Return tmp equivalent particulate scattering factor
                inline tmp<volScalarField> sigmap(const scalar dt) const;

        


        // Main calculation loop

            //- If use source correction  
//             void cellValueSourceCorrection();
            
            //- calc one cell heat transfer
            scalar calcHeatTransfer
            (
                const label celli,
                const scalar dt,
                const scalar Re,
                const scalar Pr,
                const scalar kappa,
                const scalar NCpW,
                const scalar Sh,
                const scalar Gc,
                scalar& dhsTrans,
                scalar& Sph
            );
            
            //- calc surface values for one cell 
            void calcSurfaceValues
            (
                const label celli,
                const scalar T,
                scalar& Ts,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappas
            );
            
            //- solve heat transfer for the bed
            void solveHeat(const scalar dt);
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoBedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermoBed.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
