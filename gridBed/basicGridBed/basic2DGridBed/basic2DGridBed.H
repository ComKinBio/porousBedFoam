/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::basic2DGridBed

Description
    The bed is structured, gravity is in -y direction.
    The fluid mesh is supposed to be the same as the bed grid
    (one bed grid cell -> one mesh cell)

SourceFiles
    basic2DGridBed.H
    basic2DGridBed.C

\*---------------------------------------------------------------------------*/

#ifndef basic2DGridBed_H
#define basic2DGridBed_H

#include "fvCFD.H"
#include "labelList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "DynamicList.H"
#include "demandDrivenEntry.H"
#include "typeInfo.H"
#include "ParticleForce.H"
#include "Random.H"

#include "SLGThermo.H"

// #include "ParticleForceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class integrationScheme;

template<class BedType>
class ParticleForce;

/*---------------------------------------------------------------------------*\
                         Class basic2DGridBed Declaration
\*---------------------------------------------------------------------------*/

class basic2DGridBed

{  
    
public:
        
    //- Force models type
//         typedef ParticleForce<basic2DGridBed> forceType;
        
protected:  
    
    // Protected data
    
    //initial input information
    
        //- bedName
        const word& bedName_;
    
        //fluid Mesh class
        const fvMesh& mesh_;
        
        //- Dictionary of bed properties
        IOdictionary bedProperties_;
        
        //- Dictionary of output properties
        IOdictionary outputProperties_;

        //- Dictionary of bed properties
        const dictionary particleProperties_;

        //- Dictionary of output properties
        const dictionary constProperties_;

        //- Sub-models dictionary
        const dictionary subModelProperties_;
        
    
    //- Random number generator - used by some injection routines
    mutable Random rndGen_;
        
        
    // References to the carrier gas fields
        
        //- fluid grid center list 
        const vectorField& fluidGridCellCenters_;
        
        //- fluid grid volume list 
        const scalarField& fluidGridCellVolumes_;
 
        //- gas Density [kg/m^3]
        const volScalarField& rho_;

        //- Velocity [m/s]
        const volVectorField& U_;

        //- Dynamic viscosity [Pa.s]
        const volScalarField& mu_;
        
    
    //- Gravity
        const dimensionedVector& g_;
        
    // bed submodels
        
        //- Optional particle forces
        autoPtr<ParticleForce<basic2DGridBed>> forces_;
            
        //- Momentum source term shceme
        bool momentumExplicit_;
        
        
    // bed constant properties

        //- bed min alpha fraction
        const scalar alphaMin_;
        
        //- bed material density [kg/m^3]
        demandDrivenEntry<scalar> rhop0_; 
        
        //- particle sphericity
        demandDrivenEntry<scalar> sphericity_;   
        
        
    // bed properties fields
                
        //- total particle number
        volScalarField particleNumber_;
        
        //- overall particle density
        volScalarField rhop_;

        //- volume average diameter of the particle fields [m]
        volScalarField dp_;
        
        //- surface average diameter of the particle fields [m]
        volScalarField dp2nd_;

        //- bed alpha fraction field
        volScalarField alpha_;
        
        
    // basic2DGridBed calculated members

        //- restructure the mesh, Grid[x][y] = mesh grid id, -2 means no mesh
        //- bed Grid id To Fine Grid id
        labelListList bedGridToFineGrid_;
        
        //- bed grid list, Grid[x][y] = true
        labelList bedGridList_;
        
        //- bed grid volume, read from fluid mesh 
        scalarList bedGridVolume_;
        
        //- loop listList, mesh grid list from bot to top value is fluid mesh id
        labelListList bedGridColumn_;
                
        //- bed bottom mesh list, value is fluid mesh id
        labelList bedGridBottom_;
        
        
        //- The list below is determined by the particle number
        //- should be editable within this class
        //- bed grid ID list, value is fluid mesh id
        labelList bedIDList_;
                       
        //- bed bottom mesh list, value is fluid mesh id
        labelList bedIDBottom_;
        
        //- bed top mesh list, value is fluid mesh id
        //- if the total particle number in the fluid cell > 0
        labelList bedIDTop_;
        
        //- bed internal mesh list, exclude top and bottom, value is fluid mesh id
        labelList bedIDInternal_;
        
        //- bed top&bot mesh list, value is fluid mesh id
        labelList bedIDTopBottom_;
        
        
    // Reference to the particle integration schemes

            //- Velocity integration
            autoPtr<integrationScheme> momentumIntegrator_;
        
    
    // Sources

        //- Momentum
        autoPtr<volVectorField::Internal> UTrans_;

        //- Coefficient for carrier phase U equation
        autoPtr<volScalarField::Internal> UCoeff_;
        
    
    // AddAndAbsorption
    //- (every evolve of particles)

        //- particle numbers added recored per times 
        DynamicList<scalar> addedNum_;

        //- particle numbers removed recored per times
        DynamicList<scalar> absorbedNum_;
        
        
    // Initialisation

        //- Set cloud sub-models
        void setModels();


    // bed calculate functions

        //- PreSolve 
        void preSolve();
        
        //- resetSourceTerms 
        void resetSourceTerms();


public:
        
    //- Runtime type information
    TypeName("basic2DGridBed");
  
    // Constructors

        //- Construct from fine mesh and input settings
        basic2DGridBed
        (
            const word& bedName,
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            const SLGThermo& thermo
        );
        
        //- Disallow default bitwise copy construction
        basic2DGridBed(const basic2DGridBed&) = delete;


    //- Destructor
    virtual ~basic2DGridBed()
    {}
    
    // Member Functions

        // Access
        
            // References to the mesh and databases
            
                //- Return reference to the bed name
                inline const word& bedName() const;

                //- Return reference to the mesh
                inline const fvMesh& mesh() const;

                //- Return bed properties dictionary
                inline const IOdictionary& bedProperties() const;
                
                //- Return output properties dictionary
                inline const IOdictionary& outputProperties() const;

                //- Return non-const access to the output properties dictionary
                inline IOdictionary& outputProperties();
                
                //- Return particle properties dictionary
                inline const dictionary& particleProperties() const;

                //- Return the constant properties
                inline const dictionary& constProperties() const;

                    //- Return reference to the sub-models dictionary
                inline const dictionary& subModelProperties() const;
                
            
                
            //- Return reference to the random object
            inline Random& rndGen() const;
        
                
            // References to the carrier gas fields

                //- Return carrier gas velocity
                inline const volVectorField& U() const;

                //- Return carrier gas density
                inline const volScalarField& rho() const;

                //- Return carrier gas dynamic viscosity
                inline const volScalarField& mu() const;
                
            
            
            // References to bed constant properties

                //- bed Min alpha fraction
                inline scalar alphaMin() const;
                
                //- bed material density [kg/m^3]
                inline scalar rhop0() const;
                
                //- particle sphericity
                inline scalar sphericity() const;

            
            // Environmental properties

                //- Gravi
                inline const dimensionedVector& g() const;


            //- Optional particle forces
            inline const ParticleForce<basic2DGridBed>& forces() const;

            //- Return the optional particle forces
            inline ParticleForce<basic2DGridBed>& forces();
                
                
            //- Momentum source scheme
            inline bool momentumExplicit() const;
            
        
            // Fields calculated
            
                //- Return the particle volume field [m^3]
                inline const tmp<volScalarField> bedVolume() const;
            
                //- Return the particle mass field [kg]
                inline const tmp<volScalarField> bedMass() const;

                //- Return the particle mass fraction field
                inline const tmp<volScalarField> theta() const;

                //- Return the particle volume fraction field
                inline const tmp<volScalarField> beta() const;
           
            
            // Integration schemes

                //-Return reference to velocity integration
                inline const integrationScheme& momentumIntegrator() const; 
                
            
            // Sources

                // Momentum

                    //- Return reference to momentum source
                    inline volVectorField::Internal& UTrans();

                    //- Return const reference to momentum source
                    inline const volVectorField::Internal&
                        UTrans() const;

                     //- Return coefficient for carrier phase U equation
                    inline volScalarField::Internal& UCoeff();

                    //- Return const coefficient for carrier phase U equation
                    inline const volScalarField::Internal&
                        UCoeff() const;

                    //- Return tmp momentum source term
                    inline tmp<fvVectorMatrix> SU(volVectorField& U, const scalar dt) const;
                    
            
            // Bed recored information
                
                //- Return const access to particle numbers added recored per times 
                inline scalarList addedNum() const;

                //- Return const access to particle numbers removed recored per times
                inline scalarList removedNum() const;
                
            
            // Fields protected
            
                //- Return const access to total particle number
                inline volScalarField particleNumber() const;
                
                //- Return const access to overall density
                inline volScalarField rhop() const;
                
                //- Return const access to volume average diameter of the particle fields [m]
                inline volScalarField dp() const;
                
                //- Return const access to surface average diameter of the particle fields [m]
                inline volScalarField dp2nd() const;

                //- Return const access to bed alpha fraction field
                inline volScalarField alpha() const;
                
            
            // References to the bed lists
                
                //- Return const access to bed Grid id To Fine Grid id
                inline labelListList bedGridToFineGrid() const;
                
                //- Return const access bed grid ID list
                inline labelList bedGridList() const;
                
                //- Return const access bed grid volume
                inline scalarList bedGridVolume() const;
                
                //- Return const access to bed column list
                inline labelListList bedGridColumn() const;
                        
                //- Return const access to bed bottom mesh list
                inline labelList bedGridBottom() const;
                
                
                //- Return const access bed grid ID list
                inline labelList bedIDList() const;
                
                //- Return const access to bed bottom mesh list
                //- if the total particle number in the fluid cell > 0
                inline labelList bedIDBottom() const;
                
                //- Return const access to bed top mesh list
                //- if the total particle number in the fluid cell > 0
                inline labelList bedIDTop() const;
                
                //- Return const access to bed internal mesh list
                //- exclude top and bottom, value is fluid mesh id
                inline labelList bedIDInternal() const;
                
                //- Return const access bed top&bot mesh list
                inline labelList bedIDTopBottom() const;
                    
                    
        // Edit
            
            // Fields protected

                //- Return access to total particle number
                inline volScalarField& particleNumber();
                
                //- Return access to toverall density
                inline volScalarField& rhop();
                
                //- Return access to volume average diameter of the particle fields [m]
                inline volScalarField& dp();
                
                //- Return access to surface average diameter of the particle fields [m]
                inline volScalarField& dp2nd();

                //- Return access to bed alpha fraction field
                inline volScalarField& alpha();
                
            // References to the bed lists
                
                //- updateBedList() to check current bed
                //- Return const access bed grid ID list
                inline labelList& bedIDList();
                        
                //- Return access to bed bottom mesh list
                inline labelList& bedIDBottom();
                
                //- Return access to bed top mesh list
                inline labelList& bedIDTop();
                
                //- Return access to bed internal mesh list
                //- exclude top and bottom, value is fluid mesh id
                inline labelList& bedIDInternal();
                
                //- Return access bed top&bot mesh list
                inline labelList& bedIDTopBottom();
                
                
 
                
        //- Calculate alpha field
        inline void alphaCalc();
                
                
    // Helper functions
 

        //- Particle mass
        inline scalar mass(const label celli) const;

        //- Particle volume 
        //- (TODO make particle shape as submodels in RTS)
        inline scalar volume(const label celli) const;

        //- Particle volume for a given diameter
        inline static scalar volume(const scalar dp);

        //- Particle projected area
        inline scalar areaP(const label celli) const;

        //- Projected area for given diameter
        inline static scalar areaP(const scalar dp2nd);

        //- Particle surface area
        inline scalar areaS(const label celli) const;

        //- Surface area for given diameter
        inline static scalar areaS(const scalar dp2nd);

        //- Reynolds number
        inline scalar Re(const label celli) const;

        //- Reynolds number for given conditions
        inline static scalar Re
        (
            const scalar rhoc,
            const vector& Uc,
            const scalar dp,
            const scalar muc
        );
        
        
        //- Only if bed models that included solid movement
        //- this particle of models will only be called in its inherited class through <BedType>
        //- update BedList according to current particle number field
        void updateBedList();
        
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
        
        //- Solve the bed - calculated the momentum transfer
        void solve(const scalar dt);
        
//         
//         //- update Bed overall state
//         void updateBed();      
  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#include "basic2DGridBedI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
