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
    Foam::WakaoM

Description
    The Ranz-Marshall correlation for heat transfer

\*---------------------------------------------------------------------------*/

#ifndef WakaoM_H
#define WakaoM_H

#include "MassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class WakaoM Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class WakaoM
:
    public MassTransferModel<BedType>
{

public:

    //- Runtime type information
    TypeName("Wakao");


    // Constructors

        //- Construct from dictionary
        WakaoM(const dictionary& dict, BedType& bed);

        //- Construct copy
        WakaoM(const WakaoM<BedType>& im);

        //- Construct and return a clone
        virtual autoPtr<MassTransferModel<BedType>> clone() const
        {
            return autoPtr<MassTransferModel<BedType>>
            (
                new WakaoM<BedType>(*this)
            );
        }


    //- Destructor
    virtual ~WakaoM();


    // Member Functions

        // Evaluation

            //- Sherwood number
            virtual scalar Sh
            (
                const scalar Re,
                const scalar Sc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WakaoM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
