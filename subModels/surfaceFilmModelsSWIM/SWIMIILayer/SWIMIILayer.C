/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd
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

#include "SWIMIILayer.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcFlux.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// Sub-models
#include "filmThermoModel.H"
#include "filmViscosityModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "filmRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SWIMIILayer, 0);

addToRunTimeSelectionTable(surfaceFilmRegionModel, SWIMIILayer, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

wordList SWIMIILayer::hsBoundaryTypes()
{
    wordList bTypes(T_.boundaryField().types());
    forAll(bTypes, patchi)
    {
        if
        (
            T_.boundaryField()[patchi].fixesValue()
         || isA<mixedFvPatchScalarField>(T_.boundaryField()[patchi])
         || isA<mappedFieldFvPatchField<scalar>>(T_.boundaryField()[patchi])
        )
        {
            bTypes[patchi] = fixedValueFvPatchField<scalar>::typeName;
        }
    }

    return bTypes;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool SWIMIILayer::read()
{
    // No additional properties to read
    return kinematicSingleLayer::read();
}


void SWIMIILayer::resetPrimaryRegionSourceTerms()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar(hsSp_.dimensions(), Zero);
}


void SWIMIILayer::correctThermoFields()
{
    rho_ == filmThermo_->rho();
    sigma_ == filmThermo_->sigma();
    Cp_ == filmThermo_->Cp();
    kappa_ == filmThermo_->kappa();
}


void SWIMIILayer::correctHsForMappedT()
{
    T_.correctBoundaryConditions();

    volScalarField::Boundary& hsBf = hs_.boundaryFieldRef();

    forAll(hsBf, patchi)
    {
        const fvPatchField<scalar>& Tp = T_.boundaryField()[patchi];
        if (isA<mappedFieldFvPatchField<scalar>>(Tp))
        {
            hsBf[patchi] == hs(Tp, patchi);
        }
    }
}


void SWIMIILayer::updateSurfaceTemperatures()
{
    correctHsForMappedT();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        UIndirectList<scalar>(Tw_, pp.faceCells()) =
            T_.boundaryField()[patchi];
    }
    Tw_.correctBoundaryConditions();

    // Update film surface temperature
    Ts_ = T_;
    Ts_.correctBoundaryConditions();
}


void SWIMIILayer::transferPrimaryRegionThermoFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionThermoFields();
    impSp_=rhoSp_;

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
    TParcels_.correctBoundaryConditions();
//    Info<<" min/max TParcels "<<gMin(TParcels_)<<", "<<gMax(TParcels_)<<endl;
//    Info<<" min/max Ts "<<gMin(Ts_)<<", "<<gMax(Ts_)<<endl;
//    Info<<" min/max Tw "<<gMin(Tw_)<<", "<<gMax(Tw_)<<endl;
}


void SWIMIILayer::transferPrimaryRegionSourceFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    volScalarField::Boundary& hsSpPrimaryBf =
        hsSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hsSpPrimaryBf, patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1.0/deltaT)/primaryMesh().magSf().boundaryField()[patchi]
        );

        hsSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
//    Info<<"before min, max hsSp_ "<<gMin(hsSp_)<<", "<<gMax(hsSp_)<<nl;
    hsSp_.correctBoundaryConditions();
//    Info<<"after min, max hsSp_ "<<gMin(hsSp_)<<", "<<gMax(hsSp_)<<nl;
}


void SWIMIILayer::correctAlpha()
{
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(alpha_, i)
        {
            if ((alpha_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                alpha_[i] = 1.0;
            }
            else if ((alpha_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                alpha_[i] = 0.0;
            }
        }

        alpha_.correctBoundaryConditions();
    }
    else
    {
        alpha_ ==
            pos0(delta_ - dimensionedScalar("deltaWet", dimLength, deltaWet_));
    }
}


void SWIMIILayer::updateSubmodels()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update radiation
    radiation_->correct();

    // Update injection model - mass returned is mass available for injection
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassTrans_,
        primaryEnergyTrans_
    );

    const volScalarField rMagSfDt((1/time().deltaT())/magSf());

    // Vapour recoil pressure
    pSp_ -= sqr(rMagSfDt*primaryMassTrans_)/(2*rhoPrimary_);

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, primaryMassTrans_, primaryEnergyTrans_);

    // Update source fields
    rhoSp_ += rMagSfDt*(cloudMassTrans_ + primaryMassTrans_);
    hsSp_ += rMagSfDt*(cloudMassTrans_*hs_ + primaryEnergyTrans_);
//    Info<<"min, max cloudmass "<<gMin(cloudMassTrans_)<<", "<<gMax(cloudMassTrans_)<<nl;
//    Info<<"min, max primarymass "<<gMin(primaryMassTrans_)<<", "<<gMax(primaryMassTrans_)<<nl;
//    Info<<"min, max primaryEnergy "<<gMin(primaryEnergyTrans_)<<", "<<gMax(primaryEnergyTrans_)<<nl;

    turbulence_->correct();
}


 tmp<fvScalarMatrix> SWIMIILayer::q(volScalarField& hs) const
{
    return
    (
        // Heat-transfer to the primary region
      - fvm::Sp(htcs_->h()/Cp_, hs)
      + htcs_->h()*(hs/Cp_ + alpha_*(TPrimary_ - T_))

        // Heat-transfer to the wall
      - fvm::Sp(htcw_->h()/Cp_, hs)
      + htcw_->h()*(hs/Cp_ + alpha_*(Tw_- T_))
    );
}

void SWIMIILayer::solveContinuity()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

//    Info<<"SWIM cont rhoSp_ "<<gMin(rhoSp_)<<", "<<gMax(rhoSp_)<<nl;
    solve
    (
        fvm::ddt(deltaRho_)
      + fvc::div(phi_)
     ==
      - rhoSp_
      + iceSp_
    );
}


void SWIMIILayer::solveThickness
(
    const volScalarField& pu,
    const volScalarField& pp,
    const fvVectorMatrix& UEqn
)
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    volScalarField rUA(1.0/UEqn.A());
    U_ = rUA*UEqn.H();

    surfaceScalarField deltarUAf(fvc::interpolate(delta_*rUA));
    surfaceScalarField rhof(fvc::interpolate(rho_));

    surfaceScalarField phiAdd
    (
        "phiAdd",
        regionMesh().magSf()
      * (
            fvc::snGrad(pu, "snGrad(p)")
          + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
        )
      - fvc::flux(rho_*gTan())
    );
    constrainFilmField(phiAdd, 0.0);

    surfaceScalarField phid
    (
        "phid",
        fvc::flux(U_*rho_) - deltarUAf*phiAdd*rhof
    );
    constrainFilmField(phid, 0.0);

    surfaceScalarField ddrhorUAppf
    (
        "deltaCoeff",
        fvc::interpolate(delta_)*deltarUAf*rhof*fvc::interpolate(pp)
    );

    regionMesh().setFluxRequired(delta_.name());

//    Info<<"kine thickness rhoSp_ "<<gMin(rhoSp_)<<", "<<gMax(rhoSp_)<<nl;
    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        // Film thickness equation
        fvScalarMatrix deltaEqn
        (
            fvm::ddt(rho_, delta_)
          + fvm::div(phid, delta_)
          - fvm::laplacian(ddrhorUAppf, delta_)
         ==
          - rhoSp_
	  + iceSp_
        );

        deltaEqn.solve();

        if (nonOrth == nNonOrthCorr_)
        {
            phiAdd +=
                fvc::interpolate(pp)
              * fvc::snGrad(delta_)
              * regionMesh().magSf();

            phi_ == deltaEqn.flux();
        }
    }

    // Bound film thickness by a minimum of zero
    delta_.max(0.0);

    // Update U field
    U_ -= fvc::reconstruct(deltarUAf*phiAdd);

    // Remove any patch-normal components of velocity
    U_ -= nHat()*(nHat() & U_);

    U_.correctBoundaryConditions();

    // Update film wall and surface velocities
    updateSurfaceVelocities();

    // Continuity check
    continuityCheck();
}




void SWIMIILayer::solveEnergy()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

/*
    volScalarField chi 
      = 
    volScalarField zna = hs_ + Lf_;


    iceSp_ = chi / zna;*/
/*    volScalarField imp=impSp_ * Cp_ * (TParcels_-Tf_);
    Info<<"min, max hsSp_ "<<gMin(hsSp_)<<", "<<gMax(hsSp_)<<nl;
    Info<<"min, max rhoSp_ "<<gMin(rhoSp_)<<", "<<gMax(rhoSp_)<<nl;
    Info<<"min, max hs_ "<<gMin(hs_)<<", "<<gMax(hs_)<<nl;
    Info<<"min, max zna "<<gMin(zna)<<", "<<gMax(zna)<<nl;
    Info<<"min, max imp "<<gMin(imp)<<", "<<gMax(imp)<<nl;
*/
    iceSp_=(
        - hsSp_
	+ (q(hs_) & hs_) 
        + radiation_->Shs()
	- impSp_ * Cp_ * (TParcels_-Tf_)
        +  rhoSp_ * hs_
    )
    /
    (hs_+Lf_)
    ;
    iceSp_.min(0.0);
    solve
    (
        fvm::ddt(rhoIce_, deltaIce_)
     ==
      - iceSp_
    );
//    deltaIce_.max(0.0);

    correctThermoFields();

    // Evaluate viscosity from user-model
    viscosity_->correct(pPrimary_, T_);

    // Update film wall and surface temperatures
    updateSurfaceTemperatures();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SWIMIILayer::SWIMIILayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    thermoSingleLayer(modelType, mesh, g, regionType, false),
    Tf_
    (
        IOobject
        (
            "Tf_",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("Tf_", dimTemperature, 273.16),
        zeroGradientFvPatchScalarField::typeName
    ),
    Lf_("Lf_", dimensionSet(0, 2,-2,0,0,0,0), readScalar(coeffs_.lookup("Lf"))),
    rhoIce_("rhoIce_", dimensionSet(1,-3,0,0,0,0,0), readScalar(coeffs_.lookup("rhoIce"))),
    iceSp_
    (
        IOobject
        (
            "iceSpf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimTime/dimArea, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    deltaIce_
    (
        IOobject
        (
            "deltaIcef",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    TParcels_
    (
        IOobject
        (
            "TParcels", // Same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimTemperature, Zero),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    impSp_
    (
        IOobject
        (
            "impSpf",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimTime/dimArea, Zero)
    )
{
    if (coeffs().readIfPresent("Tmin", Tmin_))
    {
        Info<< "    limiting minimum temperature to " << Tmin_ << endl;
    }

    if (coeffs().readIfPresent("Tmax", Tmax_))
    {
        Info<< "    limiting maximum temperature to " << Tmax_ << endl;
    }

    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        time().timeName(),
                        regionMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar(dimless, Zero),
                    pSp_.boundaryField().types()
                )
            );
        }
    }

    if (hydrophilic_)
    {
        coeffs_.lookup("hydrophilicDryScale") >> hydrophilicDryScale_;
        coeffs_.lookup("hydrophilicWetScale") >> hydrophilicWetScale_;
    }

    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctAlpha();

        correctThermoFields();

        // Update derived fields
        hs_ == hs(T_); //-hs(Tf_);

        deltaRho_ == delta_*rho_;

        surfaceScalarField phi0
        (
            IOobject
            (
                "phi",
                time().timeName(),
                regionMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                false
            ),
            fvc::flux(deltaRho_*U_)
        );

        phi_ == phi0;

        // Evaluate viscosity from user-model
        viscosity_->correct(pPrimary_, T_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SWIMIILayer::~SWIMIILayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void SWIMIILayer::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchi,
        facei,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    if (debug)
    {
        Info<< "    energy   = " << energySource << nl << endl;
    }

    hsSpPrimary_.boundaryFieldRef()[patchi][facei] -= energySource;
}


void SWIMIILayer::preEvolveRegion()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::preEvolveRegion();
    primaryEnergyTrans_ == dimensionedScalar(dimEnergy, Zero);
}


void SWIMIILayer::evolveRegion()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Solve continuity for deltaRho_
    solveContinuity();

    // Update sub-models to provide updated source contributions
    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=1; oCorr<=nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> tpu(this->pu());

        // Implicit pressure source coefficient
        tmp<volScalarField> tpp(this->pp());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

        // Solve energy for hs_ - also updates thermo
        solveEnergy();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn());
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update temperature using latest hs_
//    T_ == T(hs_);
}


const volScalarField& SWIMIILayer::Cp() const
{
    return Cp_;
}


const volScalarField& SWIMIILayer::kappa() const
{
    return kappa_;
}


const volScalarField& SWIMIILayer::T() const
{
    return T_;
}


const volScalarField& SWIMIILayer::Ts() const
{
    return Ts_;
}


const volScalarField& SWIMIILayer::Tw() const
{
    return Tw_;
}


const volScalarField& SWIMIILayer::hs() const
{
    return hs_;
}

const volScalarField& SWIMIILayer::deltaIce() const
{
  return deltaIce_;
}


void SWIMIILayer::info()
{
    thermoSingleLayer::info();

    Info<< indent << "min/max(deltaIce)  = "
        << gMin(deltaIce_) << ", "
        << gMax(deltaIce_) << nl
        << indent << "current ice mass   = "
        << gSum((rhoIce_*deltaIce_*magSf())()) << nl;

}


tmp<volScalarField::Internal> SWIMIILayer::Srho() const
{
    tmp<volScalarField::Internal> tSrho
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Srho",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );

    scalarField& Srho = tSrho.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchi = intCoupledPatchIDs()[i];

        scalarField patchMass =
            primaryMassTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchMass);

        const label primaryPatchi = primaryPatchIDs()[i];
        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> SWIMIILayer::Srho
(
    const label i
) const
{
    const label vapId = thermo_.carrierId(filmThermo_->name());

    tmp<volScalarField::Internal> tSrho
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Srho(" + Foam::name(i) + ")",
                time_.timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho.ref();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchi = intCoupledPatchIDs_[i];

            scalarField patchMass =
                primaryMassTrans_.boundaryField()[filmPatchi];

            toPrimary(filmPatchi, patchMass);

            const label primaryPatchi = primaryPatchIDs()[i];
            const labelUList& cells =
                primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> SWIMIILayer::Sh() const
{
    tmp<volScalarField::Internal> tSh
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Sh",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );

    scalarField& Sh = tSh.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];

        scalarField patchEnergy =
            primaryEnergyTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchEnergy);

        const label primaryPatchi = primaryPatchIDs()[i];
        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }

    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace surfaceFilmModels

// ************************************************************************* //
