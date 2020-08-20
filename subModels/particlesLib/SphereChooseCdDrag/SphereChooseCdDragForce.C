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

\*---------------------------------------------------------------------------*/

#include "SphereChooseCdDragForce.H"
//#include "Vector.H"
//#include "Tuple2.H"
//#include "volFields.H"

using namespace Foam;

template<class CloudType>
const Foam::Enum
<
    typename Foam::SphereChooseCdDragForce<CloudType>::CdType
>
Foam::SphereChooseCdDragForce<CloudType>::CdTypeNames
({
    { CdType::dfCdPutnam, "Putnam" },
    { CdType::dfCdHabashi, "Habashi" },
    { CdType::dfCdPrihodko, "Prihodko" },
    { CdType::dfCdGent, "Gent" },
    { CdType::dfCdOchkov, "Ochkov" },
    { CdType::dfCdShillerNeyman, "ShillerNeyman" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::SphereChooseCdDragForce<CloudType>::setCdType()
{
}


template<class CloudType>
Foam::scalar Foam::SphereChooseCdDragForce<CloudType>::CdRe
(
    const scalar Re
) const
{
	scalar CdRe;

	switch (CdType_)
	{
    	case CdType::dfCdPutnam:
            if (Re > 1000.0)
				CdRe = 0.424*Re;
            else
				CdRe = 24.0*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0));
            break;
        case CdType::dfCdHabashi:
            if (Re > 1300.0)
				CdRe = 0.4*Re;
            else
                CdRe = 24.0*(1.0 + 0.15*pow(Re, 0.687) + 0.00026*pow(Re, 1.38));
            break;
        case CdType::dfCdPrihodko:
            CdRe = 21.12 + 6.3*pow(Re, 0.5) + 0.25*Re;
            break;
        case CdType::dfCdGent:
            if (Re > 1000.0)
                CdRe = 0.4*Re;
            else
                CdRe = 24.0*(1.0 + 0.197*pow(Re, 0.63) + 0.00026*pow(Re, 1.38) );
            break;
        case CdType::dfCdOchkov:
            if (Re > 576.0)
                CdRe = 0.5*Re;
            else 
				if (Re > 4.0)
                	CdRe = 12.0;
            	else
					CdRe = 24.0;
            break;
        case CdType::dfCdShillerNeyman:
            if (Re > 1000.0)
                CdRe = 0.44*Re;
            else
                CdRe = 24.0*(1.0 + 0.15*pow(Re, 0.687));
            break;
        default:
        {
            FatalErrorInFunction
                << "Unhandled Cd type "
                << CdTypeNames[CdType_]
                << exit(FatalError);
        }
	}
	
//	scalar Cd;
//	Cd = CdRe/Re;
	
//  Info<< "Particle Cd = " << CdRe << endl;
//	Info<< "Particle Re = " << Re << endl;
//	Info<< "Particle Cd = " << Cd << endl;

	return CdRe;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SphereChooseCdDragForce<CloudType>::SphereChooseCdDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    CdType_(CdTypeNames.get("CdType", this->coeffs()))
{
    setCdType();
}


template<class CloudType>
Foam::SphereChooseCdDragForce<CloudType>::SphereChooseCdDragForce
(
    const SphereChooseCdDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    CdType_(df.CdType_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SphereChooseCdDragForce<CloudType>::~SphereChooseCdDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::SphereChooseCdDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    value.Sp() = mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d()));

    return value;
}


// ************************************************************************* //
