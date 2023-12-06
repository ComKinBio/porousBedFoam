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
    Foam::MassTransferModel

Description
    Templated heat transfer model class

SourceFiles
    MassTransferModel.C
    MassTransferModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef MassTransferModel_H
#define MassTransferModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "BedSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class MassTransferModel Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class MassTransferModel
:
    public BedSubModelBase<BedType>
{

public:

    //- Runtime type information
    TypeName("massTransferModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        MassTransferModel,
        dictionary,
        (
            const dictionary& dict,
            BedType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        MassTransferModel(BedType& owner);

        //- Construct from dictionary
        MassTransferModel
        (
            const dictionary& dict,
            BedType& owner,
            const word& type
        );

        //- Construct copy
        MassTransferModel(const MassTransferModel<BedType>& htm);

        //- Construct and return a clone
        virtual autoPtr<MassTransferModel<BedType>> clone() const = 0;


    //- Destructor
    virtual ~MassTransferModel();


    //- Selector
    static autoPtr<MassTransferModel<BedType>> New
    (
        const dictionary& dict,
        BedType& owner
    );


    // Member Functions

        // Evaluation

            //- Sherwood number
            virtual scalar Sh
            (
                const scalar Re,
                const scalar Sc
            ) const = 0;

            //- Return mass transfer coefficient
            virtual scalar mtc
            (
                const scalar dp,
                const scalar Re,
                const scalar Sc,
                const scalar D
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeMassTransferModel(BedType)                                         \
                                                                               \
    typedef Foam::BedType::bioBedType bioBedType;                              \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::MassTransferModel<bioBedType>,                                   \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            MassTransferModel<bioBedType>,                                     \
            dictionary                                                         \
        );                                                                     \
    }


#define makeMassTransferModelType(SS, BedType)                                 \
                                                                               \
    typedef Foam::BedType::bioBedType bioBedType;                              \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<bioBedType>, 0);              \
                                                                               \
    Foam::MassTransferModel<bioBedType>::                                      \
        adddictionaryConstructorToTable<Foam::SS<bioBedType>>                  \
            add##SS##BedType##bioBedType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MassTransferModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
