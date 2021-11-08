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
    with thermol properties and conversion

    Sub-models include:
    - particle surface heat transfer
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
    public BedType
{

public:

    // Public Typedefs

        //- Type of bed this bed was instantiated for
        typedef BedType bedType;

        //- Convenience typedef for this bed type
        typedef BioBed<BedType> BioBedType;
        
    // Static enum declear
        
        enum particlePhase 
        {
            wet_, 
            dry_, 
            char_, 
            ash_
        };
         
        static const particlePhase bedComponents[] = {wet_, dry_, char_, ash_};
    
    
    
    

    //- Class to hold BioBed constant properties
    class constantProperties
    {
        protected:

            // Protected data

                //- Constant properties dictionary
                const dictionary dict_;


        private:

            // Private Data

                //- Parcel type id - used for post-processing to flag the type
                //  of parcels issued by this cloud
                demandDrivenEntry<label> parcelTypeId_;

                //- Minimum density [kg/m^3]
                demandDrivenEntry<scalar> rhoMin_;

                //- Particle density [kg/m^3] (constant)
                demandDrivenEntry<scalar> rho0_;

                //- Minimum parcel mass [kg]
                demandDrivenEntry<scalar> minParcelMass_;


        public:

            // Constructors

                //- Null constructor
                constantProperties();

                //- Copy constructor
                constantProperties(const constantProperties& cp);

                //- Construct from dictionary
                constantProperties(const dictionary& parentDict);


            // Member Functions

                //- Return const access to the constant properties dictionary
                inline const dictionary& dict() const;

                //- Return const access to the parcel type id
                inline label parcelTypeId() const;

                //- Return const access to the minimum density
                inline scalar rhoMin() const;

                //- Return const access to the particle density
                inline scalar rho0() const;

                //- Return const access to the minimum parcel mass
                inline scalar minParcelMass() const;
    };


protected:

    // Protected data

        // Parcel properties

            //- Active flag - tracking inactive when active = false
            bool active_;

            //- Parcel type id
            label typeId_;

            //- Number of particles in Parcel
            scalar nParticle_;

            //- Diameter [m]
            scalar d_;

            //- Target diameter [m]
            scalar dTarget_;

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Density [kg/m^3]
            scalar rho_;

            //- Age [s]
            scalar age_;

            //- Time spent in turbulent eddy [s]
            scalar tTurb_;

            //- Turbulent velocity fluctuation [m/s]
            vector UTurb_;
            
    // Protected Member Functions

        //- Calculate new particle velocity
        


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

            //- Return const access to type id
            inline label typeId() const;

            //- Return const access to number of particles
            inline scalar nParticle() const;

            //- Return const access to diameter
            inline scalar d() const;

            //- Return const access to target diameter
            inline scalar dTarget() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to density
            inline scalar rho() const;

            //- Return const access to the age
            inline scalar age() const;

            //- Return const access to time spent in turbulent eddy
            inline scalar tTurb() const;

            //- Return const access to turbulent velocity fluctuation
            inline const vector& UTurb() const;


        // Edit

            //- Return const access to active flag
            inline bool& active();

            //- Return access to type id
            inline label& typeId();

            //- Return access to number of particles
            inline scalar& nParticle();

            //- Return access to diameter
            inline scalar& d();

            //- Return access to target diameter
            inline scalar& dTarget();

            //- Return access to velocity
            inline vector& U();

            //- Return access to density
            inline scalar& rho();

            //- Return access to the age
            inline scalar& age();

            //- Return access to time spent in turbulent eddy
            inline scalar& tTurb();

            //- Return access to turbulent velocity fluctuation
            inline vector& UTurb();


        // Helper functions

            //- Particle mass
            inline scalar mass() const;

            //- Particle moment of inertia around diameter axis
            inline scalar momentOfInertia() const;

            //- Particle volume
            inline scalar volume() const;

            //- Particle volume for a given diameter
            inline static scalar volume(const scalar d);

            //- Particle projected area
            inline scalar areaP() const;

            //- Projected area for given diameter
            inline static scalar areaP(const scalar d);

            //- Particle surface area
            inline scalar areaS() const;

            //- Surface area for given diameter
            inline static scalar areaS(const scalar d);


            //- Reynolds number for given conditions
            inline static scalar Re
            (
                const scalar rhoc,
                const vector& U,
                const vector& Uc,
                const scalar d,
                const scalar muc
            );

            //- Weber number for given conditions
            inline static scalar We
            (
                const scalar rhoc,
                const vector& U,
                const vector& Uc,
                const scalar d,
                const scalar sigma
            );

            //- Eotvos number for given conditions
            inline static scalar Eo
            (
                const vector& g,
                const scalar rho,
                const scalar rhoc,
                const vector& U,
                const scalar d,
                const scalar sigma
            );


        // Main calculation loop

            //- Set cell values
            
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
