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
    Foam::NoMassTransfer

Description
    Dummy heat transfer model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoMassTransfer_H
#define NoMassTransfer_H

#include "MassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NoMassTransfer Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class NoMassTransfer
:
    public MassTransferModel<BedType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoMassTransfer(const dictionary&, BedType& owner);

        //- Construct copy
        NoMassTransfer(const NoMassTransfer<BedType>& im);

        //- Construct and return a clone
        virtual autoPtr<MassTransferModel<BedType>> clone() const
        {
            return autoPtr<MassTransferModel<BedType>>
            (
                new NoMassTransfer<BedType>(*this)
            );
        }


    //- Destructor
    virtual ~NoMassTransfer();


    // Member Functions

        //- Flag to indicate whether model activates heat transfer model
        virtual bool active() const;

        //- Sherwood number
        virtual scalar Sh(const scalar, const scalar) const;

        //- Schmidt number
        virtual scalar Sc() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoMassTransfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
