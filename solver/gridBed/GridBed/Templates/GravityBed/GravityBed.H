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
    Foam::GravityBed

Description
    GravityBed class with multi-components: wet wood, dry wood, char, and ash
    with thermol properties and conversion

    Sub-models include:
    - particle surface heat transfer
    - Drying, Devo, and char reactions

SourceFiles
    GravityBedI.H
    GravityBed.C

\*---------------------------------------------------------------------------*/

#ifndef GravityBed_H
#define GravityBed_H

#include "IOstream.H"
#include "autoPtr.H"
#include "demandDrivenEntry.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class BedType>
class GravityBed;

/*---------------------------------------------------------------------------*\
                         Class GravityBed Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class GravityBed
:
    public BedType
{

public:

    // Public Typedefs

        //- Type of bed this bed was instantiated for
        typedef BedType bedType;

        //- Convenience typedef for this bed type
        typedef GravityBed<BedType> GravityBedType;
        

protected:

    // Protected data          
        
        // bed constant properties
    
            //- Dictionary of bed properties
            const dictionary collapseSetting_;

            //- bed gravity collapse max alpha fraction
            const scalar alphaGCollapse_;
            
            //- bed collapse to min alpha fraction
            const scalar alphaGCollapseMin_;
            
            //- bed collapse alpha tolerance
            const scalar alphaGCollapseTol_;



    // Protected Member Functions

        //- check whether need to calculate collapse
        bool chechCollapse();
        
        //- calculate transfer number, from j to i
        scalar collapseNumber
        (
            const label celli, 
            const label cellj
        );
        
        //- column and collapse just for once for everycell
        bool collapseOneColumn(const labelList list);

private:
    
    // Skip collapse, if method unapplicable
    bool skipCollapse_;
       

public:

    // Static Data Members

        //- Runtime type information
        TypeName("GravityBed");


    // Constructors

        //- Construct from mesh, coordinates and topology        
        //  Other properties initialised as null
        inline GravityBed
        (
            const word& bedName,
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo
        );

//         //- Construct as a copy
//         GravityBed(const GravityBed& p);
// 
//         //- Construct as a copy
//         GravityBed(const GravityBed& p, const polyMesh& mesh);

        //- Destructor
        virtual ~GravityBed()
        {}



    // Member Functions

        // Access
        
            //- Return particle properties dictionary
            inline const dictionary& collapseSetting() const;

            //- Return const access to collapse max alpha fraction
            inline const scalar& alphaGCollapse() const;

            //- Return const access to collapse to min alpha fraction
            inline const scalar& alphaGCollapseMin() const;
            
            //- Return const access to collapse alpha tolerance
            inline const scalar& alphaGCollapseTol() const;


        // Main calculation loop

            //- column and collapse 
            void gravityCollapse();
            
            //- column and collapse just for once for everycell
            void gravityCollapseOnce();
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GravityBedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GravityBed.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
