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
    Foam::RanzMarshallM

Description
    The Ranz-Marshall correlation for heat transfer

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshallM_H
#define RanzMarshallM_H

#include "MassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RanzMarshallM Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class RanzMarshallM
:
    public MassTransferModel<BedType>
{

public:

    //- Runtime type information
    TypeName("RanzMarshall");


    // Constructors

        //- Construct from dictionary
        RanzMarshallM(const dictionary& dict, BedType& bed);

        //- Construct copy
        RanzMarshallM(const RanzMarshallM<BedType>& im);

        //- Construct and return a clone
        virtual autoPtr<MassTransferModel<BedType>> clone() const
        {
            return autoPtr<MassTransferModel<BedType>>
            (
                new RanzMarshallM<BedType>(*this)
            );
        }


    //- Destructor
    virtual ~RanzMarshallM();


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
    #include "RanzMarshallM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
