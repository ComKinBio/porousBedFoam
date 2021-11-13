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
        typedef BioBed<BedType> BioBedType;
        
        using BedType::mesh_;
        
    
protected:

    // Protected data

        // References to the cloud sub-models

            // //- Mass transfer sub-model
            // autoPtr<PhaseChangeModel<ReactingCloud<CloudType>>>
            //     phaseChangeModel_;
            
            // //- Reacting phase change model
            // autoPtr<PhaseChangeModel<ReactingCloud<CloudType>>>
            //     phaseChangeModel_;

            // //- Devolatilisation model
            // autoPtr
            // <
            //     DevolatilisationModel<ReactingMultiphaseCloud<CloudType>>
            // >
            // devolatilisationModel_;

            // //- Surface reaction model
            // autoPtr
            // <
            //     SurfaceReactionModel<ReactingMultiphaseCloud<CloudType>>
            // >
            // surfaceReactionModel_;
            


        // Sources

            // //- Mass transfer fields - one per carrier phase specie
            // PtrList<volScalarField::Internal> rhoTrans_;
   
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
        TypeName("BioBed");


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline BioBed
        (
            const polyMesh& mesh,
            const vector coarseGridOrigin_,
            const vector coarseGridVertex_,
            const scalar coarseGridSize_  
        );

        //- Construct as a copy
        BioBed(const BioBed& p);

        //- Construct as a copy
        BioBed(const BioBed& p, const polyMesh& mesh);

        //- Destructor
        virtual ~BioBed();



    // Member Functions

        // Access

            //- Return const access to active flag
            inline bool active() const;

           

        // Edit

            //- Return const access to active flag
            inline bool& active();

            //- Return access to type id

        
        // Sub-models
            
            
            
        // Sources
            
            
            
        // Helper functions

            //- Particle mass
            inline scalar mass() const;

          


        // Main calculation loop

            //- Thermo conversions
            
                void calcDrying();
                
                void calcDevo();
                
                void calcCharBurnOut();
                
                //- According to the converted dmass, update bioBedFields
                //- including check exchange numbers between bio components
                void dryingUpdate();
    
                void devoUpdate(); //?  maybe one function
                
                void charBurnOutUpdate(); //?  maybe one function
                
            //- solve heat transfer and thermo conversion
            void solveThermoConversion(const scalar dt);
            
            
            //- Reload bed functions for one bio component
                
                //- Reload from ThermoBed
            
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
                
                //- Reload from basic 
                   
                    //- update Bed fields according to particle number change
                    void updateBedfieldsPrompt
                    (
                        const label fromCell,
                        const label toCell,
                        const scalar number
                    );
                    
                    //- update Bed fields according to particle number change
                    //- reload function add initial particle from source
                    void updateBedfieldsPrompt
                    (
                        const scalar number,
                        const label cell,
                        const scalar dp,
                        const scalar dp2nd
                    );
                    
                    //- update Bed fields according to particle number change
                    //- from Cell to absorption state 
                    void updateBedfieldsPrompt
                    (
                        const label cell,
                        const scalar number
                    );
                    
                    //- update Bed fields according to particle number change
                    void updateBedfields
                    (
                        const labelList& ownerCells,
                        const labelListList& ownerInteractCells,
                        const scalarListList& number // + in, - out
                    ); 
                
            
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
