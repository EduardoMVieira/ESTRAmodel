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

#include "ESTRA.H"
#include "coupledSolver.H"
#include "blockOperators.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{
    defineTypeNameAndDebug(ESTRA, 0);
    addToRunTimeSelectionTable(constitutiveEq, ESTRA, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::ESTRA::ESTRA
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    lambda_(dict.lookup("lambda")),
    L2_0_(dict.lookup("L2")),
    L2_
    (
        IOobject
        (
            "L2" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        L2_0_
    ),
    gamma_
    (
        IOobject
        (
            "gamma_" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimless, dict.lookup<scalar>("gamma"))
    ),
    tau0_
    (
        IOobject
        (
            "tau0_" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar
        (
            dimensionSet(1, -1, -2, 0, 0),
            dict.lookup<scalar>("tau0")
        )
    ),
    N1_
    (
        IOobject
        (
            "N1_" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar
        (
            dimensionSet(1, -1, -2, 0, 0),
            dict.lookup<scalar>("tau0")
        )
    ),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    A_
    (
        IOobject
        (
            "A" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor("I", L2_0_.dimensions(), symmTensor::I),
        tau_.boundaryField().types()
    ),
    varF_
    (
        IOobject
        (
            "varF_" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.0/((L2_0_ + tr(tau_)*lambda_/etaP_)/(L2_0_ - 3.0))
    ),
    iteration_(1),
    breakingInterval_(dict.lookup<label>("breakingInterval")),
    thermoLambdaPtr_(thermoFunction::New("thermoLambda", U.mesh(), dict)),
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict))
{
    checkForStab(dict);

    checkIfCoupledSolver(U.mesh().solutionDict(), tau_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constitutiveEqs::ESTRA::~ESTRA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::ESTRA::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
    // Update temperature-dependent properties
    volScalarField lambda(thermoLambdaPtr_->createField(lambda_));
    volScalarField etaP(thermoEtaPtr_->createField(etaP_));

    dimensionedSymmTensor Ist("Identity", varF_.dimensions(), symmTensor::I);

    // Velocity gradient tensor
    volTensorField L(gradU == nullptr ? fvc::grad(U())() : *gradU);

    // Convected derivative term
    volTensorField C(A_ & L);

    // First difference of normal tensions
    N1_ =
        tau_.component(tensor::XX)
      + tau_.component(tensor::YY)
      + tau_.component(tensor::ZZ);

    // Update varF
    volScalarField a(L2_/(L2_ - 3));
    gamma_ = 1.0 - mag(N1_)/tau0_;
    varF_ = 1.0 + (gamma_*(tr(A_))/L2_);

    // Stress transport equation
    fvSymmTensorMatrix AEqn
    (
        fvm::ddt(A_)
      + fvm::div(phi(), A_)
     ==
        twoSymm(C)
      - fvm::Sp(varF_/lambda, A_)
      + (a/lambda) * Ist
    );

    AEqn.relax();
    AEqn.solve();

    gamma_ = 1.0 - mag(N1_)/tau0_;
    varF_ = 1.0 + (gamma_*(tr(A_))/L2_);

    tau_ = (etaP/lambda) * (varF_*A_ - a*Ist);

    tau_.correctBoundaryConditions();

    // Breaking test
    if (std::fmod(iteration_, breakingInterval_) == 0)
    {
        forAll(A_, celli)
        {
            if ((tr(A_[celli])) > (L2_[celli]/mag(gamma_[celli])))
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
