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

Class
    Foam::SurfaceReactionModel

Description
    Templated surface reaction model class

SourceFiles
    SurfaceReactionModel.C
    SurfaceReactionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceReactionModel_H
#define SurfaceReactionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "BedSubModelBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class SurfaceReactionModel Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class SurfaceReactionModel
:
    public BedSubModelBase<BedType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //-Runtime type information
    TypeName("surfaceReactionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        SurfaceReactionModel,
        dictionary,
        (
            const dictionary& dict,
            BedType& bed
        ),
        (dict, bed)
    );


    // Constructors

        //- Construct null from owner
        SurfaceReactionModel(BedType& owner);

        //- Construct from dictionary
        SurfaceReactionModel
        (
            const dictionary& dict,
            BedType& bed,
            const word& type
        );

        //- Construct copy
        SurfaceReactionModel(const SurfaceReactionModel<BedType>& srm);

        //- Construct and return a clone
        virtual autoPtr<SurfaceReactionModel<BedType>> clone() const = 0;


    //- Destructor
    virtual ~SurfaceReactionModel();


    //- Selector
    static autoPtr<SurfaceReactionModel<BedType>> New
    (
        const dictionary& dict,
        BedType& bed
    );


    // Member Functions

        //- Update surface reactions
        //  Returns the heat of reaction
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d,
            const scalar d2nd,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar muc,
            const scalar Rec,
            const scalar e_bed,
            const scalar mass,
            const scalar massAsh,
            const scalar N,
            scalar& dmass,
            scalarField& dMassSRCarrier
        ) const = 0;

        //- Add to devolatilisation mass
        void addToSurfaceReactionMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSurfaceReactionModel(BedType)                                      \
                                                                               \
    typedef Foam::BedType::bioBedType                                          \
        bioBedType;                                                            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::SurfaceReactionModel<bioBedType>,                                \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            SurfaceReactionModel<bioBedType>,                                  \
            dictionary                                                         \
        );                                                                     \
    }


#define makeSurfaceReactionModelType(SS, BedType)                              \
                                                                               \
    typedef Foam::BedType::bioBedType                                          \
        bioBedType;                                                            \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<bioBedType>, 0);                                             \
                                                                               \
    Foam::SurfaceReactionModel<bioBedType>::                                   \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<bioBedType>>                                                 \
        add##SS##BedType##bioBedType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceReactionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
