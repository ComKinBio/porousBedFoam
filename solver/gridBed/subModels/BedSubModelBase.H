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
    Foam::BedSubModelBase

Description
    Base class for bed sub-models

SourceFiles
    BedSubModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef BedSubModelBase_H
#define BedSubModelBase_H

#include "subModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BedSubModelBase Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class BedSubModelBase
:
    public subModelBase
{
protected:

    // Protected Data

        //- Reference to the bed
        BedType& owner_;


public:

    // Constructors

        //- Construct null from owner bed
        BedSubModelBase(BedType& owner);

        //- Construct from owner bed without name
        BedSubModelBase
        (
            BedType& owner,
            const dictionary& dict,
            const word& baseName,
            const word& modelType,
            const word& dictExt = "Coeffs"
        );

        //- Construct from owner bed with name
        BedSubModelBase
        (
            const word& modelName,
            BedType& owner,
            const dictionary& dict,
            const word& baseName,
            const word& modelType
        );

        //- Copy constructor
        BedSubModelBase(const BedSubModelBase<BedType>& smb);


    //- Destructor
    virtual ~BedSubModelBase();

    //- Type of bed this model was instantiated for
    typedef BedType bedType;


    // Member Functions

        // Access

            //- Return const access to the owner bed
            const BedType& owner() const;

            //- Flag to indicate when to write a property
            virtual bool writeTime() const;


        // Edit

            //- Return non-const access to the owner bed for manipulation
            BedType& owner();


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BedSubModelBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
