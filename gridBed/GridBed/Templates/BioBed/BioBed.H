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
    Foam::BioBed

Description
    BioBed class with multi-components: wet wood, dry wood, char, and ash
    with thermol properties and conversion should be inherited from 
    thermoBed template

    Sub-models include:
    - Drying, Devo, and char reactions

SourceFiles
    BioBedI.H
    BioBed.C

\*---------------------------------------------------------------------------*/

#ifndef BioBed_H
#define BioBed_H

#include "IOstream.H"
#include "autoPtr.H"
#include "demandDrivenEntry.H"
#include "SLGThermo.H"
#include "bioBedFields.H"

// #include "ParticleForceList.H" // TODO

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class BedType>
class BioBed;


template<class bedType>
class MassTransferModel;

template<class bedType>
class DryingModel;

template<class bedType>
class DevolatilisationModel;

template<class bedType>
class SurfaceReactionModel;

/*---------------------------------------------------------------------------*\
                         Class BioBed Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class BioBed
:
    public BedType,
    public bioBedFields
{

public:

    // Public Typedefs

        //- Type of bed this bed was instantiated for
        typedef BedType bedType;

        //- Convenience typedef for this bed type
        typedef BioBed<BedType> bioBedType;
        
        using BedType::mesh_;
        
    
protected:

    // Protected data
    
         //- Constant reaction properties
    
            //- Devolatilisation activation temperature [K]
            demandDrivenEntry<scalar> TDevol_;

            //- Latent heat of devolatilisation [J/kg]
            demandDrivenEntry<scalar> LDevol_;
            
            //- shrinkageMode factor
            demandDrivenEntry<scalar> shrinkageFactorAlpha_;
            demandDrivenEntry<scalar> shrinkageFactorBeta_;
            demandDrivenEntry<scalar> shrinkageFactorGamma_;

            //- Fraction of enthalpy retained by parcel due to surface
            //  reactions
            demandDrivenEntry<scalar> hRetentionCoeff_;


        // References to the cloud sub-models

            //- Mass transfer sub-model
            autoPtr<MassTransferModel<BioBed<BedType>>>
                massTransferModel_;
            
            //- Reacting phase change model
            autoPtr<DryingModel<BioBed<BedType>>>
                dryingModel_;

            //- Devolatilisation model
            autoPtr
            <
                DevolatilisationModel<BioBed<BedType>>
            >
            devolatilisationModel_;

            //- Surface reaction model
            autoPtr
            <
                SurfaceReactionModel<BioBed<BedType>>
            >
            surfaceReactionModel_;
            
            //- Mass source term shceme
            bool massExplicit_;
            
            //- Species source term shceme
            bool speciesExplicit_;
            
            //- Double random conversion
            bool randomConvesion_;
                       

        // Sources

            //- Mass transfer fields - one per carrier phase specie
            PtrList<volScalarField::Internal> rhoTrans_;
   
   
    // Protected Member Functions

         // Initialisation

            //- Set cloud sub-models
            void setModels();
            
        
        // bed calculate functions

//             //- PreSolve 
//             void preSolve();
//             
//             //- resetSourceTerms 
//             void resetSourceTerms();


public:

    // Static Data Members

        //- Runtime type information
        TypeName("BioBed");


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        BioBed
        (
            const word& bedName,
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo
        );

//         //- Construct as a copy
//         BioBed(const BioBed& p);
// 
//         //- Construct as a copy
//         BioBed(const BioBed& p, const polyMesh& mesh);

        //- Destructor
        virtual ~BioBed()
        {}



    // Member Functions

        // Access

//             //- Return const access to active flag
//             inline bool active() const;

            //- Return const access to the devolatilisation temperature
            inline scalar TDevol() const;

            //- Return const access to the latent heat of devolatilisation
            inline scalar LDevol() const;
            
            //- Return const access to shrinkageMode factor
            inline scalar shrinkageFactorAlpha() const;
            inline scalar shrinkageFactorBeta() const;
            inline scalar shrinkageFactorGamma() const;

            //- Return const access to the fraction of enthalpy retained by
            //  parcel due to surface reactions
            inline scalar hRetentionCoeff() const;

           

        // Edit
/*
            //- Return const access to active flag
            inline bool& active();*/

            //- Return access to type id

        
        // Sub-models
            
            //- Return const reference to heat transfer model
            inline const MassTransferModel<BioBed<bedType>>&
                    massTransfer() const;
                    
            //- Return const reference to drying model
            inline const DryingModel<BioBed<bedType>>&
                    drying() const;
                    
            //- Return const reference to devolatilisation model
            inline const DevolatilisationModel<BioBed<BedType>>&
                    devolatilisation() const;
                    
            //- Return const reference to surfaceReaction model
            inline const SurfaceReactionModel<BioBed<BedType>>&
                    surfaceReaction() const;
                    
            //- Return reference to heat transfer model
            inline MassTransferModel<BioBed<bedType>>&
                    massTransfer();
                    
            //- Return reference to drying model
            inline DryingModel<BioBed<bedType>>&
                    drying();
                    
            //- Return reference to devolatilisation model
            inline DevolatilisationModel<BioBed<BedType>>&
                    devolatilisation();
                    
            //- Return reference to surfaceReaction model
            inline SurfaceReactionModel<BioBed<BedType>>&
                    surfaceReaction();
                    
            //- Mass source scheme
            inline bool massExplicit() const;
            
            //- Species source scheme
            inline bool speciesExplicit() const;
            
            //- Double random conversion
            inline bool randomConvesion() const;
            
            
        // Sources
            
            //- Mass

                //- Return reference to mass source for field i
                inline volScalarField::Internal&
                    rhoTrans(const label i);

                //- Return const access to mass source fields
                inline const PtrList<volScalarField::Internal>&
                    rhoTrans() const;

                //- Return reference to mass source fields
                inline PtrList<volScalarField::Internal>&
                    rhoTrans();

                //- Return mass source term for specie i - specie eqn
                inline tmp<fvScalarMatrix> SYi
                (
                    const label i,
                    volScalarField& Yi
                ) const;

                //- Return tmp mass source for field i - fully explicit
                inline tmp<volScalarField::Internal>
                    Srho(const label i) const;

                //- Return tmp total mass source for carrier phase
                //  - fully explicit
                inline tmp<volScalarField::Internal> Srho() const;

                //- Return total mass source term [kg/m^3/s]
                inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;
            
            
            
        // Helper functions

            //- New dp after shrinkage
            inline scalar dpAfterShrink
            (
                const scalar dp, 
                const scalar mass, 
                const scalar shrinkage
            ) const;

          


        // Main calculation loop

            //- Thermo conversions
            
                void calcDrying();
                
                void calcDevo();
                
                void calcCharBurnOut();
                
                //- Main calculation scheme functions
                //- According to the converted dmass, update bioBedFields
                //- including check exchange numbers between bio components
                void dryingUpdate();
    
                void devoUpdate(); 
                
                void charBurnOutUpdate();
                
                //- solve heat transfer for the bed
                void solveConversion(const scalar dt);
                
            //- Components update, update averaged bed properties 
            //- According to the components
            void componentUpdate();
                
            //- solve heat transfer and thermo conversion
//             void solveThermoConversion(const scalar dt);
            
            
            //- Reload bed functions for one bio component
                
                //- Reload from ThermoBed
            
                    //- calc one cell heat transfer
                    scalar calcHeatTransfer
                    (
                        const scalar dt,
                        const scalar d,
                        const scalar d2nd,
                        const scalar m,
                        const scalar T,
                        const scalar cpbed,
                        const scalar Tc,
                        const scalar Re,
                        const scalar Pr,
                        const scalar kappa,
                        const scalar NCpW,
                        const scalar Sh,
                        const scalar Gc,
                        const scalar epsilon0,
                        scalar& dhsTrans,
                        scalar& Sph
                    );
//                 
//                 //- Reload from basic 
//                    
//                     //- update Bed fields according to particle number change
//                     void updateBedfieldsPrompt
//                     (
//                         const label fromCell,
//                         const label toCell,
//                         const scalar number
//                     );
//                     
//                     //- update Bed fields according to particle number change
//                     //- reload function add initial particle from source
//                     void updateBedfieldsPrompt
//                     (
//                         const scalar number,
//                         const label cell,
//                         const scalar dp,
//                         const scalar dp2nd
//                     );
//                     
//                     //- update Bed fields according to particle number change
//                     //- from Cell to absorption state 
//                     void updateBedfieldsPrompt
//                     (
//                         const label cell,
//                         const scalar number
//                     );
//                     
//                     //- update Bed fields according to particle number change
//                     void updateBedfields
//                     (
//                         const labelList& ownerCells,
//                         const labelListList& ownerInteractCells,
//                         const scalarListList& number // + in, - out
//                     ); 
                
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BioBedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BioBed.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
