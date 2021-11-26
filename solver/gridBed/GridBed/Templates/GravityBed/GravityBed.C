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

\*---------------------------------------------------------------------------*/

#include "GravityBed.H"
#include "meshTools.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //
template<class BedType>
bool Foam::GravityBed<BedType>::chechCollapse()
{
    forAll(this->bedIDBottom_,i)
    {
        if (this->alpha_[this->bedIDBottom_[i]] > alphaGCollapse_)
        {
            return true;
        }
    }
    
    forAll(this->bedIDInternal_,i)
    {
        if (this->alpha_[this->bedIDInternal_[i]] > alphaGCollapse_)
        {
            return true;
        }
    }
    
    return false;
}


template<class BedType> Foam::scalar
Foam::GravityBed<BedType>::collapseNumber
(
    const label celli, 
    const label cellj
)
{
    scalar deltVolume = this->bedGridVolume_[celli]*(this->alpha_[celli]-alphaGCollapseMin_);
    scalar celljVp = pi/6.0*pow3(this->dp_[cellj]);
    scalar number = deltVolume/celljVp;
    number = std::round(number);
    
    if (number > this->particleNumber_[cellj])
    {
        number = this->particleNumber_[cellj];
    }
    
    scalar newAlphaCelli = 1- (this->volume(celli) + number*celljVp)/this->bedGridVolume_[celli];
    
    if (newAlphaCelli < (alphaGCollapseMin_ + alphaGCollapseTol_))
    {
        number -= 1;
    }
    
    return number;
}


template<class BedType>
void Foam::GravityBed<BedType>::collapseOneColumn(const labelList list)
{
    label listSize = list.size();

    for(label i = 0; i<listSize-1; i++)
    {
        const label celli = list[i];
        const label cellj = list[i+1];
        
        if(cellj> -2 && this->alpha_[celli]>alphaGCollapse_ 
            &&this->particleNumber_[cellj]>0)
        {
            
            scalar number = collapseNumber(celli, cellj);

            this->updateBedfieldsPrompt(cellj, celli, number);
            
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class BedType>
// Foam::GravityBed<BedType>::GravityBed
// (
//     const GravityBed<BedType>& b
// )
// :
//     BedType(b),
//     collapseSetting_(b.collapseSetting_)
//     alphaGCollapse_(b.alphaGCollapse_),
//     alphaGCollapseMin_(b.alphaGCollapseMin_)
// {}
// 
// 
// template<class BedType>
// Foam::GravityBed<BedType>::GravityBed
// (
//     const GravityBed<BedType>& b,
//     const polyMesh& mesh
// )
// :
//     BedType(b, mesh),
//     collapseSetting_(b.collapseSetting_)
//     alphaGCollapse_(b.alphaGCollapse_),
//     alphaGCollapseMin_(b.alphaGCollapseMin_)
// {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BedType>
void Foam::GravityBed<BedType>::gravityCollapse()
{
    while (chechCollapse())
    {
        const labelListList bedGridColumn = this->bedGridColumn();

        forAll(bedGridColumn, c)
        {
            const labelList columnList = bedGridColumn[c];
            
            collapseOneColumn(columnList);
        }
        
        this->updateBedList();
    }
}


template<class BedType>
void Foam::GravityBed<BedType>::gravityCollapseOnce()
{
    if (chechCollapse())
    {
        const labelListList bedGridColumn = this->bedGridColumn();
    
        forAll(bedGridColumn, c)
        {
            const labelList columnList = bedGridColumn[c];
            
            collapseOneColumn(columnList);
        }
        
        this->updateBedList();
    
    }
    else
    {
        return;
    }
}



// ************************************************************************* //
