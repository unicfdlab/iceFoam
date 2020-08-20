/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::SphereChooseCdDragForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Ergun-Wen-Yu drag model for solid spheres.

\*---------------------------------------------------------------------------*/

#ifndef SphereChooseCdDragForce_H
#define SphereChooseCdDragForce_H

#include "ParticleForce.H"
//#include "volFieldsFwd.H"
#include "Enum.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SphereChooseCdDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SphereChooseCdDragForce
:
    public ParticleForce<CloudType>
{
public:
  
		//- Cd type enumeration
    enum class CdType
    {
        dfCdPutnam,
        dfCdHabashi,
        dfCdPrihodko,
        dfCdGent,
        dfCdOchkov,
        dfCdShillerNeyman
    };

    static const Enum<CdType> CdTypeNames;

private:
    // Private Data

    // Private Member Functions

	//- Drag coefficient multiplied by Reynolds number
	scalar CdRe(const scalar Re) const;

	CdType CdType_;

	void setCdType();

public:
  
    //- Runtime type information
    TypeName("SphereChooseCdDrag");


    // Constructors

        //- Construct from mesh
        SphereChooseCdDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        SphereChooseCdDragForce(const SphereChooseCdDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new SphereChooseCdDragForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SphereChooseCdDragForce();


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SphereChooseCdDragForce.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
