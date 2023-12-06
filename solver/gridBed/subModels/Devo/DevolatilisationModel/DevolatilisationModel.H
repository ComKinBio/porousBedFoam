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
    Foam::DevolatilisationModel

Description
    Templated devolatilisation model class

SourceFiles
    DevolatilisationModel.C
    DevolatilisationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef DevolatilisationModel_H
#define DevolatilisationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "BedSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DevolatilisationModel Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class DevolatilisationModel
:
    public BedSubModelBase<BedType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //- Runtime type information
    TypeName("devolatilisationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DevolatilisationModel,
        dictionary,
        (
            const dictionary& dict,
            BedType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        DevolatilisationModel(BedType& owner);

        //- Construct from dictionary
        DevolatilisationModel
        (
            const dictionary& dict,
            BedType& owner,
            const word& type
        );

        //- Construct copy
        DevolatilisationModel(const DevolatilisationModel<BedType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DevolatilisationModel<BedType>> clone() const = 0;


    //- Destructor
    virtual ~DevolatilisationModel();


    //- Selector
    static autoPtr<DevolatilisationModel<BedType>> New
    (
        const dictionary& dict,
        BedType& owner
    );


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const scalar age,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            label& canCombust,
            scalar& dMass,
            scalarField& dMassDV,
            scalar& dMassChar
        ) const = 0;

        //- Add to devolatilisation mass
        void addToDevolatilisationMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDevolatilisationModel(BedType)                                     \
                                                                               \
    typedef Foam::BedType::bioBedType                                          \
        bioBedType;                                                            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::DevolatilisationModel<bioBedType>,                               \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            DevolatilisationModel<bioBedType>,                                 \
            dictionary                                                         \
        );                                                                     \
    }


#define makeDevolatilisationModelType(SS, BedType)                             \
                                                                               \
    typedef Foam::BedType::bioBedType                                          \
        bioBedType;                                                            \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<bioBedType>, 0);                                             \
                                                                               \
    Foam::DevolatilisationModel<bioBedType>::                                  \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<bioBedType>>                                                 \
        add##SS##BedType##bioBedType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DevolatilisationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
