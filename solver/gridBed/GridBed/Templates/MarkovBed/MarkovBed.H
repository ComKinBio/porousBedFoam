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
    Foam::MarkovBed

Description
    MarkovBed class with multi-components: wet wood, dry wood, char, and ash
    with thermol properties and conversion

    Sub-models include:
    - particle surface heat transfer
    - Drying, Devo, and char reactions

SourceFiles
    MarkovBedI.H
    MarkovBed.C

\*---------------------------------------------------------------------------*/

#ifndef MarkovBed_H
#define MarkovBed_H

#include "IOstream.H"
#include "autoPtr.H"
#include "demandDrivenEntry.H"

// #include "ParticleForceList.H" // TODO

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class BedType>
class MarkovBed;

/*---------------------------------------------------------------------------*\
                         Class MarkovBed Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class MarkovBed
:
    public BedType
{

public:

    // Public Typedefs

        //- Type of bed this bed was instantiated for
        typedef BedType bedType;

        //- Convenience typedef for this bed type
        typedef MarkovBed<BedType> MarkovBedType;

    


protected:

    // Protected data

        // Parcel properties

            //- Active flag - tracking inactive when active = false
            bool active_;

            //- Parcel type id
            label typeId_;

            //- Number of particles in Parcel
            scalar nParticle_;


    // Protected Member Functions

        //- Calculate new particle velocity
        


public:

    // Static Data Members

        //- Runtime type information
        TypeName("MarkovBed");


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline MarkovBed
        (
            const polyMesh& mesh,
            const vector coarseGridOrigin_,
            const vector coarseGridVertex_,
            const scalar coarseGridSize_  
        );

        //- Construct as a copy
        MarkovBed(const MarkovBed& p);

        //- Construct as a copy
        MarkovBed(const MarkovBed& p, const polyMesh& mesh);

        //- Destructor
        virtual ~MarkovBed();



    // Member Functions

        // Access

            //- Return const access to active flag
            inline bool active() const;

            //- Return const access to type id
            inline label typeId() const;

            //- Return const access to number of particles
            inline scalar nParticle() const;


        // Edit

            //- Return const access to active flag
            inline bool& active();

            //- Return access to type id
            inline label& typeId();

            //- Return access to number of particles
            inline scalar& nParticle();


        // Helper functions


            //- Particle volume
            inline scalar volume() const;

            //- Particle volume for a given diameter
            inline static scalar volume(const scalar d);


        // Main calculation loop

            //- Set cell values
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MarkovBedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MarkovBed.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
