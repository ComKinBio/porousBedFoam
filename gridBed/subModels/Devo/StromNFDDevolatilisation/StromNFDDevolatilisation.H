/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::StromNFDDevolatilisation

Description
    Di Blasi devolatisation model.
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = Ax.exp(- Ex/R.T)
		        gas 	kappa1
        wood = 	tar	    kappa2
		        char	kappa3
		
	tar = gas	kappa4
	      char	kappa5

    Where:
        kappa = rate constant
        A1    = activation energy (user input)
        E     = pre-exponential factor (user input)
        R     = universal gas constant
        T     = temperature

    Usage:

        StromNFDDevolatilisationCoeffs
        {
            devolKinetic1	() // (A1 E1)
            devolKinetic2	() // (A2 E2)
            devolKinetic3	() // (A3 E3)
            devolKinetic4	() // (A4 E4)
            devolKinetic5	() // (A5 E5)
            volgas
            (
                (CH4     0.5 )   // (name weight_fraction)
                (CO2     0.5)   // (name weight_fraction)
            );
            
            volatileResidualCoeff 1e-6;
        }

\*---------------------------------------------------------------------------*/

#ifndef StromNFDDevolatilisation_H
#define StromNFDDevolatilisation_H

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class StromNFDDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class BedType>
class StromNFDDevolatilisation
:
    public DevolatilisationModel<BedType>
{
    // Helper class to store specie-local volatile data
    class volgas
    {
        // Private data

            //- Specie name
            word name_;

            //- Weight fraction
            scalar Y_;

    public:

        // Constructors

            //- Null constructor
            volgas()
            :
                name_(word::null),
                Y_(0.0)
            {}

            //- Construct from Istream
            volgas(Istream& is)
            :
                name_(is),
                Y_(readScalar(is))
            {}

        //- Destructor
        ~volgas()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the weight fraction
                scalar Y() const
                {
                    return Y_;
                }


        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volgas& vd)
            {
                is.readBeginList("volgas");
                is  >>  vd.name_  >> vd.Y_;
                is.readEndList("volgas");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volgas& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.Y_ 
                    << token::END_LIST;

                return os;
            }
    };

    // Private data

        // Model constants
    
            //- List of volgas1 data - (name weight fraction for volatile gas 1)
            List<volgas> volgas_;

            //- Mapping between local and bed volatile gaseous species 1
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;
            
            //- Kinetic parameters for the 1st devolatilisation reaction (A E)
            const List<scalar> devolKinetic1_;
            
            //- Kinetic parameters for the 2nd devolatilisation reaction (A E)
            const List<scalar> devolKinetic2_;
            
            //- Kinetic parameters for the 3rd devolatilisation reaction (A E)
            const List<scalar> devolKinetic3_;

            //- label for tar, should be named as "tar"
            const label id_tar;
	 

public:

    //- Runtime type information
    TypeName("StromNFDDevolatilisation");


    // Constructors

        //- Construct from dictionary
        StromNFDDevolatilisation
        (
            const dictionary& dict,
            BedType& owner
        );

        //- Construct and return a clone
        virtual autoPtr<DevolatilisationModel<BedType> > clone() const
        {
            return autoPtr<DevolatilisationModel<BedType> >
            (
                new StromNFDDevolatilisation<BedType>(*this)
            );
        }


    //- Destructor
    virtual ~StromNFDDevolatilisation();


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
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "StromNFDDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
