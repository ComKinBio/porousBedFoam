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
    Foam::HeatTransferModel

Description
    Templated heat transfer model class

SourceFiles
    HeatTransferModel.C
    HeatTransferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef HeatTransferModel_H
#define HeatTransferModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "BedSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class HeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class HeatTransferModel
:
    public BedSubModelBase<BedType>
{
    // Private Data

        //- Apply Bird's correction to the htc
        const Switch BirdCorrection_;


public:

    //- Runtime type information
    TypeName("heatTransferModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HeatTransferModel,
        dictionary,
        (
            const dictionary& dict,
            BedType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        HeatTransferModel(BedType& owner);

        //- Construct from dictionary
        HeatTransferModel
        (
            const dictionary& dict,
            BedType& owner,
            const word& type
        );

        //- Construct copy
        HeatTransferModel(const HeatTransferModel<BedType>& htm);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<BedType>> clone() const = 0;


    //- Destructor
    virtual ~HeatTransferModel();


    //- Selector
    static autoPtr<HeatTransferModel<BedType>> New
    (
        const dictionary& dict,
        BedType& owner
    );


    // Member Functions

        // Access

            //- Return the Bird htc correction flag
            const Switch& BirdCorrection() const;


        // Evaluation

            //- Nusselt number
            virtual scalar Nu
            (
                const scalar Re,
                const scalar Pr
            ) const = 0;

            //- Return heat transfer coefficient
            virtual scalar htc
            (
                const scalar dp,
                const scalar Re,
                const scalar Pr,
                const scalar kappa,
                const scalar NCpW
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeHeatTransferModel(BedType)                                       \
                                                                               \
    typedef Foam::BedType::thermoBedType thermoBedType;                  \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::HeatTransferModel<thermoBedType>,                              \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            HeatTransferModel<thermoBedType>,                                \
            dictionary                                                         \
        );                                                                     \
    }


#define makeHeatTransferModelType(SS, BedType)                               \
                                                                               \
    typedef Foam::BedType::thermoBedType thermoBedType;                  \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<thermoBedType>, 0);         \
                                                                               \
    Foam::HeatTransferModel<thermoBedType>::                                 \
        adddictionaryConstructorToTable<Foam::SS<thermoBedType>>             \
            add##SS##BedType##thermoBedType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HeatTransferModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
