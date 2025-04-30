/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025 UFES
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

#include "ESTRALog.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{
    defineTypeNameAndDebug(ESTRALog, 0);
    addToRunTimeSelectionTable(constitutiveEq, ESTRALog, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::ESTRALog::ESTRALog
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    ESTRA(name, U, phi, dict),
    theta_
    (
        IOobject
        (
            "theta" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    eigVals_
    (
        IOobject
        (
            "eigVals" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor("I", dimless, pTraits<tensor>::I),
        extrapolatedCalculatedFvPatchField<tensor>::typeName
    ),
    eigVecs_
    (
        IOobject
        (
            "eigVecs" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor("I", dimless, pTraits<tensor>::I),
        extrapolatedCalculatedFvPatchField<tensor>::typeName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constitutiveEqs::ESTRALog::~ESTRALog()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::ESTRALog::correct
(
    const volScalarField* alpha,
    const volTensorField* gradU
)
{
    // Update temperature-dependent properties
    volScalarField lambda = thermoLambdaPtr_->createField(lambda_);
    volScalarField etaP = thermoEtaPtr_->createField(etaP_);

    // Decompose grad(U).T()
    volTensorField L(gradU == nullptr ? fvc::grad(U())() : *gradU);

    dimensionedScalar c1("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0);

    volTensorField B = c1*eigVecs_;
    volTensorField omega = B;
    volTensorField M = (eigVecs_.T() & L.T() & eigVecs_);

    decomposeGradU(M, eigVals_, eigVecs_, omega, B);

    // Solve the constitutive equation in theta = log(c)

    dimensionedTensor Itensor
    (
        "Identity",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        tensor::I
    );

    N1_ =
        tau_.component(tensor::XX)
      + tau_.component(tensor::YY)
      + tau_.component(tensor::ZZ);

    scalar a = L2_0_.value()/(L2_0_.value() - 3.0);
    gamma_ = 1.0 - mag(N1_)/tau0_;

    volScalarField f =
        1.0 + (gamma_*(tr(eigVecs_ & eigVals_ & eigVecs_.T()))/L2_);

    fvSymmTensorMatrix thetaEqn
    (
        fvm::ddt(theta_)
      + fvm::div(phi(), theta_)
     ==
        symm
        (
            (omega & theta_)
          - (theta_ & omega)
          + 2.0*B
          + (1.0/lambda)
           *(
                eigVecs_ & (a*inv(eigVals_) - f*Itensor) & eigVecs_.T()
            )
        )
    );

    thetaEqn.relax();
    thetaEqn.solve();

    // Diagonalization of theta
    calcEig(theta_, eigVals_, eigVecs_);

    // Recalculate gamma
    gamma_ = 1.0 - mag(N1_)/tau0_;
    f = 1.0 + (gamma_*(tr(eigVecs_ & eigVals_ & eigVecs_.T()))/L2_);

    // Convert from theta to tau
    tau_ =
        (etaP/lambda)*symm(f*(eigVecs_ & eigVals_ & eigVecs_.T()) - a*Itensor);

    tau_.correctBoundaryConditions();

    // Breaking test
    if (std::fmod(iteration_, breakingInterval_) == 0)
    {
        forAll(theta_, celli)
        {
            if
            (
                (tr(eigVecs_[celli] & eigVals_[celli] & eigVecs_[celli].T()))
              > (L2_[celli]/mag(gamma_[celli]))
            )
            {
                L2_[celli] = (L2_[celli])*1.0404;
            }
        }
    }

    // L2 equation for advection
    fvScalarMatrix L2Eqn
    (
        fvm::ddt(L2_)
     ==
      - fvm::div(phi(), L2_)
    );

    L2Eqn.relax();
    L2Eqn.solve();

    iteration_++;
}


// ************************************************************************* //
