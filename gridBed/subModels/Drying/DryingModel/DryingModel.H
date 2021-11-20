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
    Foam::DryingModel

Description
    Templated phase change model class

SourceFiles
    DryingModel.C
    DryingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef DryingModel_H
#define DryingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "BedSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class DryingModel Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class DryingModel
:
    public BedSubModelBase<BedType>
{
public:

    // Public enumerations

        //- Enthalpy transfer type
        enum enthalpyTransferType
        {
            etLatentHeat,
            etEnthalpyDifference
        };

        //- Name representations of enthalpy transfer types
        static const Foam::wordList enthalpyTransferTypeNames;


protected:

    // Protected data

        //- Enthalpy transfer type enumeration
        enthalpyTransferType enthalpyTransfer_;


        // Counters

            //- Mass of lagrangian phase converted
            scalar dMass_;


    // Protected Member Functions

        //- Convert word to enthalpy transfer type
        enthalpyTransferType wordToEnthalpyTransfer(const word& etName) const;

public:

    //- Runtime type information
    TypeName("dryingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DryingModel,
        dictionary,
        (
            const dictionary& dict,
            BedType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        DryingModel(BedType& owner);

        //- Construct from dictionary
        DryingModel
        (
            const dictionary& dict,
            BedType& owner,
            const word& type
        );

        //- Construct copy
        DryingModel(const DryingModel<BedType>& pcm);

        //- Construct and return a clone
        virtual autoPtr<DryingModel<BedType>> clone() const = 0;


    //- Destructor
    virtual ~DryingModel();


    //- Selector
    static autoPtr<DryingModel<BedType>> New
    (
        const dictionary& dict,
        BedType& owner
    );


    // Access

        //- Return the enthalpy transfer type enumeration
        const enthalpyTransferType& enthalpyTransfer() const;


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const label celli,
            const scalar Re,
            const scalar Nu,
            const scalar Sh,
            const scalar d,
            const scalar nu,
            const scalar T,
            const scalar Ts,
            const scalar pc,
            const scalar Tc,
            const scalarField& X,
            scalarField& dMassPC
        ) const = 0;

        //- Return the enthalpy per unit mass
        virtual scalar dh
        (
            const label idc,
            const label idl,
            const scalar p,
            const scalar T
        ) const;

        //- Return vapourisation temperature
        virtual scalar Tvap(const scalarField& X) const;

        //- Return maximum/limiting temperature
        virtual scalar TMax(const scalar p, const scalarField& X) const;

        //- Add to phase change mass
        void addToDryingMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDryingModel(BedType)                                               \
                                                                               \
    typedef Foam::BedType::bioBedType bioBedType;                              \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::DryingModel<bioBedType>,                                         \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            DryingModel<bioBedType>,                                           \
            dictionary                                                         \
        );                                                                     \
    }


#define makeDryingModelType(SS, BedType)                                       \
                                                                               \
    typedef Foam::BedType::bioBedType bioBedType;                              \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<bioBedType>, 0);              \
                                                                               \
    Foam::DryingModel<bioBedType>::                                            \
        adddictionaryConstructorToTable<Foam::SS<bioBedType>>                  \
            add##SS##BedType##bioBedType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DryingModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
