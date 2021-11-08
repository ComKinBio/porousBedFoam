/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
// 
\*---------------------------------------------------------------------------*/

#include "basic2DGridBed.H"
#include "integrationScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basic2DGridBed, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::basic2DGridBed::basic2DGridBed
(
    const word& bedName,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    const SLGThermo& thermo
)
:       
    bedName_(bedName),
    mesh_(mesh),
    bedProperties_
    (
        IOobject
        (
            bedName + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    outputProperties_
    (
        IOobject
        (
            bedName + "OutputProperties",
            mesh_.time().timeName(),
            "bedOutputProperties"/bedName,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    particleProperties_(bedProperties_.subOrEmptyDict("particleProperties")),
    constProperties_(bedProperties_.subOrEmptyDict("constantProperties")),
    subModelProperties_(bedProperties_.subOrEmptyDict("subModels")),
    fluidGridCellCenters_(mesh.cellCentres()),
    fluidGridCellVolumes_(mesh.cellVolumes()),
    rho_(rho),
    U_(U),
    mu_(mu),
    g_(g),
    forces_(nullptr),
    momentumExplicit_(true),
    alphaMin_(constProperties_.lookupOrDefault<scalar>("alphaMin", 0.01)),
    rhop_(particleProperties_, "rhop"),
    particleNumber_
    (
        IOobject
        (
            IOobject::groupName(bedName, "particleNumber"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphac",dimless, 0)
    ),
    dp_
    (
        IOobject
        (
            IOobject::groupName(bedName, "dp"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, 0)
    ),
    dp2nd_
    (
        IOobject
        (
            IOobject::groupName(bedName, "dp2nd"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, 0)
    ),
    alpha_
    (
        IOobject
        (
            IOobject::groupName(bedName, "alpha"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphac",dimless, 1.0)
    ),
    bedGridToFineGrid_(0),
    bedGridList_(0),
    bedGridVolume_(0),
    bedGridColumn_(0),
    bedGridBottom_(0),
    bedIDList_(0),
    bedIDBottom_(0),
    bedIDTop_(0),
    bedIDInternal_(0),
    bedIDTopBottom_(0),
    momentumIntegrator_(nullptr),
    UTrans_
    (
        new volVectorField::Internal
        (
            IOobject
            (
                bedName + ":UTrans",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector(dimMass*dimVelocity, Zero)
        )
    ),
    UCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                bedName + ":UCoeff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMass, 0)
        )
    ),
    addedNum_(0),
    absorbedNum_(0)
{
    //- read scheme
    Istream& is = subModelProperties_.lookup("MomentumSource");
    const word scheme(is);
    if (scheme == "semiImplicit")
    {
         momentumExplicit_ = false;
    }
    else if (scheme == "explicit")
    {
        momentumExplicit_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid scheme " << scheme << ". Valid schemes are "
            << "explicit and semiImplicit" << exit(FatalError);
    }

    
    Info << "Construct coarse grid"<< endl;
        
    const dictionary bedSettings = bedProperties_.subDict("bedSettings");
    
    //Two point to decide bed grid geometry, the origin point
    const vector bedGridOrigin_ = bedSettings.lookup("gridOrigin");

    //The second point to determine bed grid geometry
    const vector bedGridVertex_ = bedSettings.lookup("gridVertex");

    //bed grid size   
    const scalar bedGridSize_ = readScalar(bedSettings.lookup("gridSize"));

    //- bed grid number in x direction
    const label bedGridNumberInX = floor(mag(bedGridVertex_.x()-bedGridOrigin_.x())/bedGridSize_);
        
    //- bed grid number in y direction
    const label bedGridNumberInY = floor(mag(bedGridVertex_.y()-bedGridOrigin_.y())/bedGridSize_);

    //iteration index
    label i,j;
    
    //creat bed grid
//     label bedGridNumber_ = bedGridNumberInX*bedGridNumberInY*bedGridNumberInZ;
    
    DynamicList<labelList> bedGridIDList2DTem(0);
    DynamicList<labelList> bedGridColumnTem(0);
    DynamicList<label> bedGridIDListTem(0);
    DynamicList<scalar> bedGridVolumeTem(0);
    DynamicList<label> bedIDBottomTem(0);

    
    for (i=1;i<=bedGridNumberInX;i++) //x loop
    {
        
        // For G[x][y]
        DynamicList<label> bedGrid2DIDListInOneColumnTem(0);
                
        // Get labelListList bedGridColumn and labelList bedGridBottom
        // bedIDBottom bedIDTOP...updated from them according to bed fields
        DynamicList<label> bedGridIDListInOneColumnTem(0);
        
        bool findCellFirstTimeInColume = true;
                
        for (j=1;j<=bedGridNumberInY;j++) //column loop
        {
            
            // check whether the cell center in the bed grid cell
            scalar xl_ = bedGridOrigin_.x() + bedGridSize_*(i-1);
            scalar xr_ = bedGridOrigin_.x() + bedGridSize_*(i);
            scalar yl_ = bedGridOrigin_.y() + bedGridSize_*(j-1); 
            scalar yr_ = bedGridOrigin_.y() + bedGridSize_*(j);
            
            bool findCell = false;
            
            //TODO remove one each time will be more efficient?
            forAll(fluidGridCellCenters_, idF) 
            {
                scalar fluidCellX_ = fluidGridCellCenters_[idF].x();
                scalar fluidCellY_ = fluidGridCellCenters_[idF].y();
     
                if
                ( 
                    (xl_ <= fluidCellX_) && (fluidCellX_ < xr_) 
                 && (yl_ <= fluidCellY_) && (fluidCellY_ < yr_)
                 && !findCell
                )
                {
                    findCell = true;
                    bedGridIDListTem.append(idF);
                    bedGridVolumeTem.append(fluidGridCellVolumes_[idF]);
                    bedGridIDListInOneColumnTem.append(idF);
                    bedGrid2DIDListInOneColumnTem.append(idF);
                    
                    if (findCellFirstTimeInColume)
                    {
                        bedIDBottomTem.append(idF);
                        findCellFirstTimeInColume = false;
                    }
                }
                
            }//all fluid grid loop
            
            if (!findCell)
            {
                bedGrid2DIDListInOneColumnTem.append(-2);
            }
            
        }//column loop finish
        
        // update [x]=[i] for labelListList
        bedGridIDList2DTem.append(bedGrid2DIDListInOneColumnTem);
        
        if (!findCellFirstTimeInColume)
        {
            bedGridColumnTem.append(bedGridIDListInOneColumnTem);
        }
        
        
    }//x loop finish
    
    
    bedGridToFineGrid_ = bedGridIDList2DTem;
    bedGridList_ = bedGridIDListTem;
    bedGridVolume_ = bedGridVolumeTem;
    bedGridColumn_ = bedGridColumnTem;
    bedGridBottom_ = bedIDBottomTem;
    
    //- update bedList
    updateBedList();
    
    setModels();
    
    alphaCalc();
    
    Info << "Bed model: "<<bedName_ <<" Constructed"<< endl;

}
    
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::basic2DGridBed::setModels()
{
    dictionary forcesModel_ = subModelProperties_.subOrEmptyDict("particleForces");
    
    forAllConstIter(IDLList<entry>, forcesModel_, iter)
    {
        const word& model = iter().keyword();
        if (iter().isDict())
        {
            forces_.reset
            (
                ParticleForce<basic2DGridBed>::New
                (
                    *this,
                    mesh_,
                    iter().dict(),
                    model
                ).ptr()
            );
        }
        else
        {
            forces_.reset
            (
                ParticleForce<basic2DGridBed>::New
                (
                    *this,
                    mesh_,
                    dictionary::null,
                    model
                ).ptr()
            );
        }
    }
    
    dictionary integrationSchemeDict_ = subModelProperties_.subOrEmptyDict("integrationSchemes");
    
    momentumIntegrator_.reset
    (
        integrationScheme::New
        (
            "Momentum",
            integrationSchemeDict_
        ).ptr()
    );
}


void Foam::basic2DGridBed::solve(const scalar dt)
{
    preSolve();
    
    //- Loop bed calculate drag force
    forAll(bedIDList_, i)
    {
        const label celli = bedIDList_[i];
        const scalar npi = particleNumber_[celli];
        const scalar massi = mass(celli);
        vector Udp = Zero;
        
        // Reynolds number and gas phase
        const scalar Rei = Re(celli);
        const scalar muc = mu_[celli];
        const vector Uc = U_[celli];

        
        // Sources
        //~~~~~~~~

        // Linearised momentum source coefficient
        scalar Spu = 0.0;

        // Momentum transfer from the particle to the carrier phase
        vector dUTrans = Zero;
        
        const forceSuSp Fcp = forces().calcCoupled(dt, massi, Rei, muc);
        
        const vector acp = (Fcp.Sp()*Uc + Fcp.Su())/massi;
        const scalar bcp = Fcp.Sp()/massi;
        
        const vector deltaU = momentumIntegrator().delta(Udp, dt, acp, bcp);

        dUTrans -= massi*deltaU;
        
        Spu = dt*Fcp.Sp();
        
        // Update momentum transfer
        UTrans()[celli] = npi*dUTrans;

        // Update momentum transfer coefficient
        UCoeff()[celli] += npi*Spu;
    }
    
}



void Foam::basic2DGridBed::preSolve()
{
    resetSourceTerms();
}
      
      
void Foam::basic2DGridBed::resetSourceTerms()
{
    UTrans().field() = Zero;
    UCoeff().field() = 0.0;
}


void Foam::basic2DGridBed::updateBedList()
{
    DynamicList<label> bedIDListTem(0);
    DynamicList<label> bedIDBottomListTem(0);
    DynamicList<label> bedIDTopListTem(0);
    DynamicList<label> bedIDInternalListTem(0);
    DynamicList<label> bedIDTopBottomListTem(0);
    
    forAll(bedGridColumn_, clmi)
    {
        labelList ithColumn = bedGridColumn_[clmi];
        label columnH = ithColumn.size()-1;
        bool isBottom = true;
        
        forAll(ithColumn, j)
        {
            scalar number = particleNumber_[ithColumn[j]];
            
            if (j < columnH)
            {
                scalar numberUp = particleNumber_[ithColumn[j+1]];
                
                if (isBottom && numberUp < 1 && number >= 1)
                {
                    bedIDListTem.append(ithColumn[j]);
                    bedIDTopBottomListTem.append(ithColumn[j]);
                    isBottom = false;
                }
                else if (isBottom && numberUp >= 1 && number >= 1)
                {
                    bedIDListTem.append(ithColumn[j]);
                    bedIDBottomListTem.append(ithColumn[j]);
                    isBottom = false;
                }
                else if (!isBottom && numberUp >= 1 && number >= 1)
                {
                    bedIDListTem.append(ithColumn[j]);
                    bedIDInternalListTem.append(ithColumn[j]);
                    isBottom = false;
                }
                else if (!isBottom && numberUp < 1 && number >= 1)
                {
                    bedIDListTem.append(ithColumn[j]);
                    bedIDTopListTem.append(ithColumn[j]);
                    isBottom = false;
                }
            }
            else if (isBottom && number >= 1)
            {
                bedIDListTem.append(ithColumn[j]);
                bedIDTopBottomListTem.append(ithColumn[j]);
                isBottom = false;
            }
            else if (!isBottom && number >= 1)
            {
                bedIDListTem.append(ithColumn[j]);
                bedIDTopListTem.append(ithColumn[j]);
                isBottom = false;
            }
        }
    }

    bedIDList_ = bedIDListTem;
    bedIDBottom_ = bedIDBottomListTem;
    bedIDTop_ = bedIDTopListTem;
    bedIDInternal_ = bedIDInternalListTem;
    bedIDTopBottom_ = bedIDTopBottomListTem;
    
}


void Foam::basic2DGridBed::updateBedfieldsPrompt
(
    const label fromCell,
    const label toCell,
    const scalar number
)
{
    // Before transfer number records
    scalar npj = particleNumber_[toCell];
    
    // update number
    particleNumber_[fromCell] -= number;
    particleNumber_[toCell] += number;
    
    // update dp_, toCell only
    dp_[toCell] = cbrt((pow3(dp_[toCell])*npj + pow3(dp_[fromCell])*number)\
                            /particleNumber_[toCell]);
                
    // update dp2nd_, toCell only
    dp2nd_[toCell] = sqrt((sqr(dp2nd_[toCell])*npj + sqr(dp2nd_[fromCell])*number)\
                            /particleNumber_[toCell]);
    
    // update alpha in i,j
    alpha_[fromCell] = max(1.0 - volume(fromCell)/mesh_.V()[fromCell], alphaMin_);
    alpha_[toCell] = max(1.0 - volume(toCell)/mesh_.V()[toCell], alphaMin_);
    
}



//- add initial particle from source
void Foam::basic2DGridBed::updateBedfieldsPrompt
(
    const scalar number,
    const label cell,
    const scalar dp,
    const scalar dp2nd
)
{
     // update dp_
    dp_[cell] = cbrt((pow3(dp_[cell])*particleNumber_[cell] + pow3(dp)*number)\
                            /(particleNumber_[cell] + number));
    
    // update dp2nd_
    dp2nd_[cell] = sqrt((sqr(dp2nd_[cell])*particleNumber_[cell] + sqr(dp2nd)*number)\
                            /(particleNumber_[cell] + number));
    
    // update number
    particleNumber_[cell] += number;
    
    // update alpha
    alpha_[cell] = max(1.0 - volume(cell)/mesh_.V()[cell], alphaMin_);
    
    addedNum_.append(number);
}



//- from Cell to absorption state
void Foam::basic2DGridBed::updateBedfieldsPrompt
(
    const label cell,
    const scalar number
)
{
    // update number
    particleNumber_[cell] -= number;
    
    // update alpha
    alpha_[cell] = max(1.0 - volume(cell)/mesh_.V()[cell], alphaMin_);
    
    absorbedNum_.append(number);
}


void Foam::basic2DGridBed::updateBedfields
(
    const labelList& ownerCells,
    const labelListList& ownerInteractCells,
    const scalarListList& number // + in, - out
)
{
    //- Save old fields
    volScalarField dpSaved = dp_;
    volScalarField dp2ndSaved = dp2nd_;

    forAll(ownerCells, i)
    {
        label celli = ownerCells[i];
        labelList iInterac = ownerInteractCells[i];
        scalarList numList = number[i];
        
        forAll(iInterac, j)
        {
            label cellj = iInterac[j];
            scalar numij = numList[j];
            
            if (numij > 0)
            {
                dp_[celli] = cbrt((pow3(dp_[celli])*particleNumber_[celli] + pow3(dpSaved[cellj])*numij)\
                            /(particleNumber_[celli] + numij));
                dp2nd_[celli] = sqrt((sqr(dp_[celli])*particleNumber_[celli] + sqr(dp2ndSaved[cellj])*numij)\
                            /(particleNumber_[celli] + numij));
                particleNumber_[celli] += numij;
            }
            else
            {
                dp_[celli] = cbrt((pow3(dp_[celli])*particleNumber_[celli] + pow3(dpSaved[celli])*numij)\
                            /(particleNumber_[celli] + numij));
                dp2nd_[celli] = sqrt((sqr(dp_[celli])*particleNumber_[celli] + sqr(dp2ndSaved[celli])*numij)\
                            /(particleNumber_[celli] + numij));
                particleNumber_[celli] += numij;
            }
            
        }
            
    }
    
    //- bed alpha fraction field
    alphaCalc();
}


// ************************************************************************* //
